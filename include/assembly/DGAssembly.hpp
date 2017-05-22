//
// Created by tyler on 5/2/17.
//

#ifndef YAFEL_DGASSEMBLY_HPP
#define YAFEL_DGASSEMBLY_HPP

#include "yafel_globals.hpp"
#include "element/ElementFactory.hpp"
#include "fe_system/FESystem.hpp"
#include "assembly/AssemblyRequirement.hpp"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <vector>


YAFEL_NAMESPACE_OPEN


/**
 * Assembly procedure for Discontinuous Galerkin FEM.
 * The procedure involves looping over mesh faces and
 * computing fluxes between elements.
 *
 * Designed for use with explicit time stepping,
 *
 * @tparam Physics
 * @param feSystem
 * @param requirements
 */
template<typename Physics>
void DGAssembly(FESystem &feSystem,
                Physics &physics,
                std::vector<AssemblyRequirement> requirements = {AssemblyRequirement::Residual})
{

    constexpr int simulation_dimension = Physics::nsd();

    //Unpack the FESystem
    auto &GlobalResidual = feSystem.getGlobalResidual();
    auto &GlobalSolution = feSystem.getSolution();
    GlobalResidual *= 0;
    auto &dofm = feSystem.getDoFManager();
    auto dof_per_node = dofm.dof_per_node;
    double time = feSystem.currentTime();

    bool assemble_tangent{false};
    bool assemble_residual{false};
    bool assemble_dt_mass{false};
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
            case AssemblyRequirement::DtDtMass:
                assemble_dtdt_mass = true;
                break;
        }
    }

    if (physics.mass_constructed)
        assemble_dt_mass = false;
    else if (assemble_dt_mass && !physics.mass_constructed) {
        physics.inverse_mass_matrices.resize(dofm.nCells());
        physics.mass_constructed = true;
    }


//#pragma omp parallel shared(GlobalResidual, tangent_triplets, dofm)
    {
        //Define thread-local variables

        //storage buffers
        std::vector<double> local_tangent_buffer;
        std::vector<double> local_dtmass_buffer;
        std::vector<double> local_residual_buffer;
        std::vector<double> local_solution_buffer;
        std::vector<int> global_dof_buffer_l;
        std::vector<int> global_dof_buffer_r;
        std::vector<int> left_local_nodes;
        std::vector<int> right_local_nodes;
        std::vector<int> face_nodes;
        std::vector<Eigen::Triplet<double>> local_triplets;
        Eigen::VectorXd private_residual = Eigen::VectorXd::Constant(GlobalResidual.rows(), 0.0);

        // Need two different element factories because
        // only a single Element of a given type is instantiated
        // in either one. In order to do element boundary flux
        // calculations, you need two separate elements simultaneously.
        ElementFactory EF_L(dofm.dof_per_node);
        ElementFactory EF_R(dofm.dof_per_node);


        // Boundary Fluxes
//#pragma omp single
        {
            for (int fi = 0; fi < dofm.interior_faces.size(); ++fi) {
                auto &F = dofm.interior_faces[fi];
                if (F.left < 0 || F.right < 0) {
                    int el{-1};
                    int fl{-1};
                    int fr{-1};
                    int rl_idx{-1};
                    if (F.left >= 0) {
                        el = F.left;
                        //fl = F.left_flocal;
                        //fr = F.left_rot;
                        //rl_idx = 0;
                        dofm.getLeftFaceNodes(fi,face_nodes);
                    } else {
                        el = F.right;
                        //fl = F.right_flocal;
                        //fr = F.right_rot;
                        //rl_idx = 1;
                        dofm.getRightFaceNodes(fi,face_nodes);
                    }

                    auto &E = EF_L.getElement(dofm.element_types[el]);
                    dofm.getGlobalDofs(el, global_dof_buffer_l);
                    //auto &face_nodes = E.face_perm[fl][fr][rl_idx];

                    for (int fqpi = 0; fqpi < E.nFQP(); ++fqpi) {
                        //auto nl = E.face_update<Physics::nsd()>(el, fqpi, F, dofm);
                        auto nl = E.face_update<Physics::nsd()>(el, fqpi, face_nodes, dofm);
                        coordinate<> xqp;
                        double U{0};
                        for (int i = 0; i < face_nodes.size(); ++i) {
                            xqp += dofm.dof_nodes[global_dof_buffer_l[face_nodes[i]]] * E.boundaryShapeValues[fqpi](i);
                            U += GlobalSolution(global_dof_buffer_l[face_nodes[i]]) * E.boundaryShapeValues[fqpi](i);
                        }
                        if (rl_idx == 1) {
                            //nl = -nl;
                        }

                        double fluxVal = Physics::BoundaryFlux(nl, xqp, time, U);

                        for (int A = 0; A < face_nodes.size(); ++A) {
                            double val = E.boundaryShapeValues[fqpi](A) * fluxVal * E.jxw;
                            GlobalResidual(global_dof_buffer_l[face_nodes[A]]) -= val;
                        }

                    }

                } else {

                    int e_left = F.left;
                    int fl_l = F.left_flocal;
                    //int fr_l = F.left_rot;

                    int e_right = F.right;
                    int fl_r = F.right_flocal;
                    //int fr_r = F.right_rot;

                    dofm.getGlobalDofs(e_left, global_dof_buffer_l);
                    dofm.getGlobalDofs(e_right, global_dof_buffer_r);

                    auto &EL = EF_L.getElement(dofm.element_types[e_left]);
                    auto &ER = EF_R.getElement(dofm.element_types[e_right]);

                    //auto &left_nodes = EL.face_perm[fl_l][fr_l][0];
                    //auto &right_nodes = EL.face_perm[fl_r][fr_r][1];

                    dofm.getLeftFaceNodes(fi,left_local_nodes);
                    dofm.getRightFaceNodes(fi,right_local_nodes);

                    Eigen::VectorXd R_l = Eigen::VectorXd::Constant(left_local_nodes.size(), 0.0);
                    Eigen::VectorXd R_r = Eigen::VectorXd::Constant(right_local_nodes.size(), 0.0);

                    for (int fqpi = 0; fqpi < EL.nFQP(); ++fqpi) {
                        //auto nl = EL.face_update<Physics::nsd()>(e_left, fqpi, F, dofm);
                        //auto nr = ER.face_update<Physics::nsd()>(e_right, fqpi, F, dofm);
                        auto nl = EL.face_update<Physics::nsd()>(e_left, fqpi, left_local_nodes, dofm);
                        auto nr = ER.face_update<Physics::nsd()>(e_right, fqpi, right_local_nodes, dofm);

                        coordinate<> xqp{0, 0, 0};
                        double Uleft{0};
                        double Uright{0};

                        for (int i = 0; i < left_local_nodes.size(); ++i) {
                            double shapeVal = EL.boundaryShapeValues[fqpi](i);
                            xqp += dofm.dof_nodes[global_dof_buffer_l[left_local_nodes[i]]] * shapeVal;//EL.boundaryShapeValues[fqpi](i);
                            Uleft += GlobalSolution(global_dof_buffer_l[left_local_nodes[i]]) * shapeVal;
                            Uright += GlobalSolution(global_dof_buffer_r[right_local_nodes[i]]) * shapeVal;
                        }

                        double fluxVal = Physics::Flux(nl, xqp, time, Uleft, Uright);

                        for (int i = 0; i < left_local_nodes.size(); ++i) {
                            int li = left_local_nodes[i];
                            int gli = global_dof_buffer_l[li];
                            int ri = right_local_nodes[i];
                            int gri = global_dof_buffer_r[ri];

                            auto shapeVal = EL.boundaryShapeValues[fqpi](i);
                            GlobalResidual(gli) -= shapeVal * fluxVal * EL.jxw;
                            GlobalResidual(gri) += shapeVal * fluxVal * ER.jxw;
                        }
                    }


                }
            }
        }


        // Element-level fluxes
//#pragma omp for
        for (int elnum = 0; elnum < dofm.nCells(); ++elnum) {

            auto et = dofm.element_types[elnum];
            if (et.topoDim != simulation_dimension) {
                continue;
            }
            auto &E = EF_L.getElement(et);
            dofm.getGlobalDofs(elnum, global_dof_buffer_l);

            auto local_dofs = E.localMesh.nNodes() * dof_per_node;
            if (assemble_tangent && static_cast<int>(local_tangent_buffer.size()) < local_dofs * local_dofs) {
                local_tangent_buffer.resize(local_dofs * local_dofs, 0.0);
            }
            if (assemble_dt_mass && static_cast<int>(local_dtmass_buffer.size()) < local_dofs * local_dofs) {
                local_dtmass_buffer.resize(local_dofs * local_dofs, 0.0);
            }
            if (assemble_residual && static_cast<int>(local_residual_buffer.size()) < local_dofs) {
                local_residual_buffer.resize(local_dofs, 0.0);
            }
            if (static_cast<int>(local_solution_buffer.size()) < local_dofs) {
                local_solution_buffer.resize(local_dofs, 0.0);
            }

            for (auto &x : local_tangent_buffer) {
                x = 0;
            }
            for (auto &x : local_dtmass_buffer) {
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
            //Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> local_dt_mass(
            //local_dtmass_buffer.data(), local_dofs, local_dofs);

            Eigen::MatrixXd local_dt_mass;
            if(assemble_dt_mass) {
                local_dt_mass = Eigen::MatrixXd::Constant(local_dofs, local_dofs, 0.0);
            }

            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> local_residual(local_residual_buffer.data(),
                                                                                local_dofs);
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> local_solution(local_solution_buffer.data(),
                                                                                local_dofs);

            //copy element solution values into local_solution
            for (int i = 0; i < global_dof_buffer_l.size(); ++i) {
                local_solution(i) = GlobalSolution(global_dof_buffer_l[i]);
                local_residual(i) = GlobalResidual(global_dof_buffer_l[i]);
            }


            auto nqp = E.nQP();
            for (int qpi = 0; qpi < nqp; ++qpi) {

                coordinate<> xqp;
                for(int A=0; A<global_dof_buffer_l.size(); ++A) {
                    xqp += dofm.dof_nodes[global_dof_buffer_l[A]]*E.shapeValues[qpi](A);
                }

                E.update<Physics::nsd()>(elnum, qpi, dofm);

                if (assemble_tangent) {
                    //Physics::LocalTangent(E, qpi, time, local_tangent);
                }
                if (assemble_residual) {
                    Physics::LocalResidual(E, qpi, xqp, time, local_solution, local_residual);
                }
                if (assemble_dt_mass) {
                    Physics::LocalMass(E, qpi, time, local_dt_mass);
                }
            }

            if(assemble_dt_mass){
                Eigen::PartialPivLU<Eigen::MatrixXd> MLU(local_dt_mass);
                //physics.inverse_mass_matrices.push_back(MLU);
                physics.inverse_mass_matrices[elnum] = MLU;
            }


            //Invert local mass matrix and solve the local_residual
            auto &MLU = physics.inverse_mass_matrices[elnum];
            local_residual = MLU.solve(local_residual);

            //Assemble into global
            for (auto A : IRange(0, local_dofs)) {
                auto GA = global_dof_buffer_l[A];
                if (assemble_tangent) {
                    for (auto B : IRange(0, local_dofs)) {
                        auto GB = global_dof_buffer_l[B];
                        local_triplets.emplace_back(GA, GB, local_tangent(A, B));
                    }
                }
                if (assemble_residual) {
                    GlobalResidual(GA) += local_residual(A);
                }
            }

        } //end elnum loop



    } // end parallel block


}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_DGASSEMBLY_HPP
