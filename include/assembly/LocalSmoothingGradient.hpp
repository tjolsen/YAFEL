//
// Created by tyler on 6/12/17.
//

#ifndef YAFEL_LOCALSMOOTHINGGRADIENT_HPP
#define YAFEL_LOCALSMOOTHINGGRADIENT_HPP

#include "yafel_globals.hpp"

#include "element/Element.hpp"
#include "element/ElementFactory.hpp"
#include "fe_system/FESystem.hpp"
#include "assembly/AssemblyRequirement.hpp"

#include <eigen3/Eigen/Core>
#include <vector>

YAFEL_NAMESPACE_OPEN

template<typename Physics>
void LocalSmoothingGradient(FESystem &feSystem)
{
    constexpr int NSD = Physics::nsd();
    constexpr int simulation_dimension = NSD;
    auto &Solution = feSystem.getSolution();
    auto &SolutionGradient = feSystem.getSolutionGradient();
    auto &dofm = feSystem.getDoFManager();

    SolutionGradient.resize(Solution.rows(), NSD);

    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> VolTimesGrad =
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::Constant(Solution.rows(), NSD, 0.0);
    Eigen::VectorXd Volume = Eigen::VectorXd::Constant(Solution.rows(), 0.0);

    ElementFactory EF;
    std::vector<int> global_dof_buffer;
    std::vector<double> local_solution_buffer;
    Eigen::MatrixXd qp_grad = Eigen::MatrixXd::Constant(dofm.dof_per_node, NSD, 0.0);

    for (int elnum = 0; elnum < dofm.nCells(); ++elnum) {

        auto et = dofm.element_types[elnum];
        if (et.topoDim != simulation_dimension) {
            continue;
        }

        auto &E = EF.getElement(et);
        dofm.getGlobalDofs(elnum, global_dof_buffer);

        auto n_local_dofs = E.localMesh.nNodes() * dofm.dof_per_node;
        if (local_solution_buffer.size() < n_local_dofs) {
            local_solution_buffer.resize(n_local_dofs, 0);
        }

        Eigen::Map<Eigen::VectorXd> local_solution(local_solution_buffer.data(), n_local_dofs);
        for(int i=0; i<n_local_dofs; ++i) {
            local_solution(i) = Solution(global_dof_buffer[i]);
        }

        auto nqp = E.nQP();

        for (auto qpi : IRange(0, nqp)) {
            qp_grad *= 0;
            E.update<NSD>(elnum, qpi, dofm);

            for (auto A: IRange(0, n_local_dofs)) {
                auto comp = E.getComp(A);
                auto node = E.getNode(A);

                for (auto i : IRange(0, NSD)) {
                    qp_grad(comp, i) += local_solution(A) * E.shapeGrad(node, i);
                }
            }

            for(auto A : IRange(0,n_local_dofs)) {
                Volume(global_dof_buffer[A]) += E.shapeValues[qpi](E.getNode(A))*E.jxw;

                for(auto i : IRange(0,NSD)) {
                    VolTimesGrad(global_dof_buffer[A], i) += E.shapeValues[qpi](E.getNode(A)) * qp_grad(E.getComp(A),i) * E.jxw;
                }
            }
        }
    }


    for(auto A : IRange(0,static_cast<int>(VolTimesGrad.rows()))) {
        for(auto i : IRange(0,NSD)) {
            auto vol = Volume(A);
            auto val = VolTimesGrad(A,i)/Volume(A);
            SolutionGradient(A,i) = val;
        }
    }

}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_LOCALSMOOTHINGGRADIENT_HPP
