//
// Created by tyler on 4/16/17.
//

#ifndef YAFEL_CGASSEMBLY_HPP
#define YAFEL_CGASSEMBLY_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "element/ElementFactory.hpp"
#include "fe_system/FESystem.hpp"
#include "assembly/AssemblyRequirement.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

YAFEL_NAMESPACE_OPEN

/**
 *
 * \brief General-purpose Continuous-Galerkin finite element assembly
 *
 * @tparam Physics Class that defines the local element matrix/vector construction in static void methods
 */
template<typename Physics>
void CGAssembly(FESystem &feSystem,
                std::vector<AssemblyRequirement> requirements = {AssemblyRequirement::Residual,
                                                                 AssemblyRequirement::Tangent})
{

    // Unpack the FESystem
    auto &GlobalTangent = feSystem.getGlobalTangent();
    auto &GlobalResidual = feSystem.getGlobalResidual();
    auto &GlobalSolution = feSystem.getSolution();
    auto &dofm = feSystem.getDoFManager();
    auto dof_per_node = dofm.dof_per_node;
    constexpr int simulation_dimension = Physics::nsd();
    auto time = feSystem.currentTime();


    bool assemble_tangent{false};
    bool assemble_residual{false};
    bool assemble_dt_mass{false};
    bool inverse_lumped_dt_mass{false};
    bool assemble_dtdt_mass{false};
    for (auto req : requirements) {
        switch (req) {
            case AssemblyRequirement::Residual:
                assemble_residual = true;
                break;
            case AssemblyRequirement::Tangent:
                assemble_tangent = true;
                break;
            case AssemblyRequirement::DtMass:
                assemble_dt_mass = true;
                break;
            case AssemblyRequirement::LumpedDtMassInverse:
                assemble_dt_mass = true;
                inverse_lumped_dt_mass = true;
                break;
            case AssemblyRequirement::DtDtMass:
                assemble_dtdt_mass = true;
                break;
        }
    }

    //shared objects
    std::vector<Eigen::Triplet<double>> tangent_triplets;
    int total_triplets{0};

#pragma omp parallel shared(tangent_triplets, GlobalResidual, dofm)
    {// open parallel block
        //storage buffers
        std::vector<double> local_tangent_buffer;
        std::vector<double> local_residual_buffer;
        std::vector<double> local_solution_buffer;
        std::vector<int> global_dof_buffer;
        Eigen::VectorXd private_residual = Eigen::VectorXd::Constant(GlobalResidual.rows(), 0.0);
        std::vector<Eigen::Triplet<double>> local_triplets;

        //Create an ElementFactory
        ElementFactory EF;
#pragma omp for
        for (int elnum = 0; elnum < dofm.nCells(); ++elnum) {

            auto et = dofm.element_types[elnum];
            if (et.topoDim != simulation_dimension) {
                continue;
            }

            auto &E = EF.getElement(et);
            dofm.getGlobalDofs(elnum, global_dof_buffer);


            auto local_dofs = E.localMesh.nNodes() * dof_per_node;
            if (assemble_tangent && static_cast<int>(local_tangent_buffer.size()) < local_dofs * local_dofs) {
                local_tangent_buffer.resize(local_dofs * local_dofs, 0.0);
            }
            if (assemble_residual && static_cast<int>(local_residual_buffer.size()) < local_dofs) {
                local_residual_buffer.resize(local_dofs, 0.0);
                local_solution_buffer.resize(local_dofs, 0.0);
            }

            for (auto &x : local_tangent_buffer) {
                x = 0;
            }
            for (auto &x : local_residual_buffer) {
                x = 0;
            }
            for (auto &x : local_solution_buffer) {
                x = 0;
            }

            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> local_tangent(
                    local_tangent_buffer.data(), local_dofs, local_dofs);
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> local_residual(local_residual_buffer.data(),
                                                                                local_dofs);
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> local_solution(local_solution_buffer.data(),
                                                                                local_dofs);

            //Fill the local solution buffer
            for(int i=0; i<global_dof_buffer.size(); ++i) {
                local_solution(i) = GlobalSolution(global_dof_buffer[i]);
            }

            auto nqp = E.nQP();
            for (auto qpi : IRange(0, nqp)) {
                coordinate<> xqp;
                for(int A=0; A<global_dof_buffer.size(); ++A) {
                    xqp += dofm.dof_nodes[global_dof_buffer[A]]*E.shapeValues[qpi](A);
                }

                E.update<Physics::nsd()>(elnum, qpi, dofm);

                if (assemble_tangent) {
                    Physics::LocalTangent(E, qpi, xqp, time, local_solution, local_tangent);
                }
                if (assemble_residual) {
                    Physics::LocalResidual(E, qpi, xqp, time, local_solution, local_residual);
                }

            }//end quadrature point loop


            //Assemble into global
            for (auto A : IRange(0, local_dofs)) {
                auto GA = global_dof_buffer[A];
                if (assemble_tangent || assemble_dt_mass || assemble_dtdt_mass) {
                    for (auto B : IRange(0, local_dofs)) {
                        auto GB = global_dof_buffer[B];
                        local_triplets.emplace_back(GA, GB, local_tangent(A, B));
                    }
                }
                if (assemble_residual) {
                    private_residual(GA) += local_residual(A);
                }
            }

        }// end element loop


#pragma omp critical
        {
            total_triplets += static_cast<int>(local_triplets.size());
        }//end critical block
#pragma omp barrier
#pragma omp single
        {
            tangent_triplets.reserve(total_triplets);
        }
#pragma omp barrier
#pragma omp critical
        {
            for (auto &trip : local_triplets) {
                tangent_triplets.push_back(trip);
            }

            GlobalResidual += private_residual;
        }
    }//end parallel block

    GlobalTangent.setFromTriplets(tangent_triplets.begin(), tangent_triplets.end());
}

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_CGASSEMBLY_HPP
