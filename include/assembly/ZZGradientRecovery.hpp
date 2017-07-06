//
// Created by tyler on 7/6/17.
//

#ifndef YAFEL_ZZGRADIENTRECOVERY_HPP
#define YAFEL_ZZGRADIENTRECOVERY_HPP

#include "yafel_globals.hpp"
#include "element/ElementFactory.hpp"
#include "fe_system/FESystem.hpp"
#include "utils/DoFManager.hpp"

#include <eigen3/Eigen/QR>
#include <vector>
#include <set>

YAFEL_NAMESPACE_OPEN

/**
 * \brief Implementation of the Zienkiwicz & Zhu gradient recovery technique
 *
 * Implementation of "ZZ" patch recovery technique for 2D and 3D linear elements.
 * It uses the quadrature points as the sampling points rather than the "superconvergent"
 * points suggested in the paper, but initial experiments indicate that it works very well.
 * In the future, the Elements class may be augmented with a "RecoveryPoints" list, to be
 * used for this procedure instead of the quadrature points. For simplex elements, this is the
 * quadrature rule that you would use for a (P-1)-order element.
 *
 * The method consists of building a node-element adjacency list for each node in the mesh,
 * then fitting a local polynomial of order P (currently, P=1) at the node using the values
 * at the quadrature points in the adjacent elements.
 */
template<typename Physics>
void ZZGradientRecovery(FESystem &fesystem)
{

    constexpr int NSD = Physics::nsd();
    auto &dofm = fesystem.getDoFManager();
    auto &Solution = fesystem.getSolution();
    auto &SolutionGradient = fesystem.getSolutionGradient();
    SolutionGradient.resize(Solution.rows(), NSD);

    //Temporary assertion. Need to generalize
    assert(dofm.polyOrder == 1 && "ZZ Recovery currently only intended for P=1 elements");

    ElementFactory EF;

    std::vector<std::set<int>> adjacent_elements(dofm.nNodes());

    // Build node-element patch adjacency
    std::vector<int> node_container;
    std::vector<int> dof_container;

    for (int elnum = 0; elnum < dofm.nCells(); ++elnum) {
        dofm.getGlobalNodes(elnum, node_container);
        for (auto n : node_container) {
            adjacent_elements[n].insert(elnum);
        }
    }



    // Loop over nodes
    for (int nodeNum = 0; nodeNum < dofm.nNodes(); ++nodeNum) {

        auto const &adj_elems = adjacent_elements[nodeNum];
        int num_sample_points{0};
        int num_basis_funcs{0};
        for (auto e : adj_elems) {
            auto et = dofm.element_types[e];
            if (et.topoDim != NSD) {
                continue;
            }
            num_sample_points += EF.getElement(et).nQP();
        }

        if (NSD == 2) {
            num_basis_funcs = 3;
        } else if (NSD == 3) {
            num_basis_funcs = 4;
        }

        //Build system: P*a = b, where P contains polynomial entries for a least-squares fit of patch values
        // EG: P = [1 x y] for 2d linear elements
        Eigen::MatrixXd P(num_sample_points, num_basis_funcs);
        Eigen::MatrixXd b(num_sample_points, NSD);
        int idx{0};
        for (auto e : adj_elems) {

            //auto dofpn = dofm.dof_per_node;
            auto et = dofm.element_types[e];
            if (et.topoDim != NSD) {
                continue;
            }
            auto &E = EF.getElement(et);
            dofm.getGlobalNodes(e, node_container);
            dofm.getGlobalDofs(e, dof_container);

            for (int qpi = 0; qpi < E.nQP(); ++qpi) {

                coordinate<> xqp;
                E.update<NSD>(e, qpi, dofm);
                Tensor<NSD,1,double> field_grad(0);


                //compute point location and field gradient
                for (int A = 0; A < E.localMesh.nNodes(); ++A) {
                    xqp += E.shapeValues[qpi](A) * dofm.dof_nodes[node_container[A]];
                    field_grad += Solution(dof_container[A])*make_TensorMap<NSD,1>(&E.shapeGrad(A,0));
                }

                P(idx, 0) = 1;
                for (int i = 0; i < NSD; ++i) {
                    P(idx, i + 1) = xqp(i);
                    b(idx,i) = field_grad(i);
                }

                ++idx;
            }
        }//end adj_elems loop


        Eigen::MatrixXd a = P.colPivHouseholderQr().solve(b);;

        // Recover gradF(x_n) = P(x_n)*a
        auto x_n = dofm.dof_nodes[nodeNum];
        Eigen::MatrixXd Pn(1,num_basis_funcs);
        Pn(0,0) = 1;
        for(int i=0; i<NSD; ++i) {
            Pn(0,i+1) = x_n(i);
        }

        Eigen::MatrixXd recovered_grad = Pn*a;

        for(int i=0; i<NSD; ++i) {
            SolutionGradient(nodeNum,i) = recovered_grad(0,i);
        }

    }

}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_ZZGRADIENTRECOVERY_HPP
