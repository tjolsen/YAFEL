//
// Created by tyler on 10/12/17.
//

#ifndef YAFEL_SUPG_ADVECTION_HPP
#define YAFEL_SUPG_ADVECTION_HPP

#include "yafel.hpp"

YAFEL_NAMESPACE_OPEN

template<int NSD>
class AdvectionDiffusionSUPG : PDEBase<NSD>
{
public:
    static constexpr int nsd() { return NSD; }

    static void LocalResidual(const Element &E, int qpi, double,
                              Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &R_el)
    {
        //Sign convention because moving it over to other side of equation
        auto dof_per_elem = E.dofPerNode()*E.localMesh.nNodes();

        coordinate<> xqp(0);
        for(int i=0; i<E.localMesh.nNodes(); ++i) {
            xqp += E.shapeValues[qpi](i)*E.localMesh.getGeometryNodes()[i];
        }
        auto U = advectionVelocity(xqp);
        U = (alpha*h/2)*U/norm(U);


        for(int A=0; A<dof_per_elem; ++A) {
            auto Anode = E.getNode(A);
            auto Acomp = E.getComp(A);

            R_el(A) -= (E.shapeValues[qpi](Anode)*Q + U(Acomp)*E.shapeGrad(Anode,Acomp))*E.jxw;
        }

    }

    static void LocalTangent(const Element &E, int qpi, double,
                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &K_el)
    {

        auto dof_per_elem = E.dofPerNode()*E.localMesh.nNodes();

        coordinate<> xqp(0);
        for(int i=0; i<E.localMesh.nNodes(); ++i) {
            xqp += E.shapeValues[qpi](i)*E.localMesh.getGeometryNodes()[i];
        }


        //advection velocity at quadrature point
        auto U = advectionVelocity(xqp);

        //Total diffusion tensor
        Tensor<NSD,2> K = TensorEye<NSD,double>()*k + (alpha*h/2)*otimes(U,U)/norm(U);

        for(int A=0; A<dof_per_elem; ++A) {
            auto Anode = E.getNode(A);
            auto Acomp = E.getComp(A);

            for(int B=0; B<dof_per_elem; ++B) {
                auto Bnode = E.getNode(B);
                auto Bcomp = E.getComp(B);

                K_el(A,B) += (
                        E.shapeValues[qpi](Anode)*U(Bcomp)*E.shapeGrad(Bnode,Bcomp) //advection term
                        + E.shapeGrad(Anode,Acomp)*E.shapeGrad(Bnode,Bcomp)*K(Acomp,Bcomp) //diffusion term (includes supg)
                )*E.jxw;

            }
        }

    }


    constexpr static double Q{0.0};
    constexpr static double k{0.01};
    constexpr static double alpha = 1.0;
    constexpr static double h = 0.1;

    static Tensor<NSD,1,double> advectionVelocity(coordinate<> x) {
        Tensor<NSD,1,double> ret(0); ret(0) = 1;
        return ret;
    };
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_SUPG_ADVECTION_HPP
