//
// Created by tyler on 10/12/17.
//

#ifndef YAFEL_SUPG_ADVECTION_HPP
#define YAFEL_SUPG_ADVECTION_HPP

#include "yafel.hpp"
#include <cmath>

YAFEL_NAMESPACE_OPEN

enum TransientFlag
{
    Steady,
    Transient
};

template<int NSD, TransientFlag flag = Steady>
class AdvectionDiffusionSUPG : PDEBase<NSD>
{
public:
    static constexpr int nsd() { return NSD; }

    template<typename TU, typename TR>
    static void LocalResidual(const Element &E, int qpi, coordinate<> const &xqp, double,
                              Eigen::DenseBase<TU> &U_el, //Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &U_el,
                              Eigen::DenseBase<TR> &R_el) //Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &R_el)
    {
        //Sign convention because moving it over to other side of equation

        Tensor<NSD, 1> gradPhi(0);
        double alpha = alphaImpl();

        if constexpr (flag == Transient)
        {
            for (int Anode = 0; Anode < E.localMesh.nNodes(); ++Anode) {
                for (int i = 0; i < NSD; ++i) {
                    gradPhi(i) += E.shapeGrad(Anode, i) * U_el(Anode);
                }
            }
        }
        //Actual velocity
        auto U = advectionVelocity(xqp);

        //Total diffusion tensor
        Tensor<NSD, 2> K = TensorEye<NSD, double>() * k + (alpha * h / 2) * otimes(U, U) / norm(U);

        //Diffusion residual
        Tensor<NSD, 1> diffusionR = K * gradPhi;

        //Rescale velocity
        auto normU = norm(U);

        Tensor<NSD, 1> Uhat = (normU>0) ? ((alpha * h / 2) * U / normU).eval() : (U*0).eval();

        for (int Anode = 0; Anode < E.localMesh.nNodes(); ++Anode) {
            for (int Acomp = 0; Acomp < NSD; ++Acomp) {

                R_el(Anode) -= (E.shapeValues[qpi](Anode) * Q + Uhat(Acomp) * E.shapeGrad(Anode, Acomp)) * E.jxw;

                if constexpr (flag == Transient)
                {
                    R_el(Anode) -= (E.shapeValues[qpi](Anode) * U(Acomp) * gradPhi(Acomp) +
                                E.shapeGrad(Anode, Acomp) * diffusionR(Acomp)) * E.jxw;
                }
            }
        }

    }

    template<typename TU>
    static void LocalTangent(const Element &E, int qpi, coordinate<> const &xqp, double,
                             Eigen::DenseBase<TU> &U_el, //Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &U_el,
                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &K_el)
    {

        //advection velocity at quadrature point
        auto U = advectionVelocity(xqp);
        double alpha = alphaImpl();

        auto normU = norm(U);

        //Total diffusion tensor
        Tensor<NSD, 2> K = TensorEye<NSD, double>() * k;
        if(normU > 0) {
            K += (alpha * h / (normU*2)) * otimes(U, U) ;
        }

        for (int Anode = 0; Anode < E.localMesh.nNodes(); ++Anode) {
            for(int Acomp=0; Acomp < NSD; ++Acomp) {

                for (int Bnode = 0; Bnode < E.localMesh.nNodes(); ++Bnode) {
                    for(int Bcomp=0; Bcomp < NSD; ++Bcomp) {

                        K_el(Anode, Bnode) += (
                                              E.shapeValues[qpi](Anode) * U(Bcomp) *
                                              E.shapeGrad(Bnode, Bcomp) //advection term
                                              + E.shapeGrad(Anode, Acomp) * E.shapeGrad(Bnode, Bcomp) *
                                                K(Acomp, Bcomp) //diffusion term (includes supg)
                                      ) * E.jxw;
                    }
                }
            }
        }

    }


    constexpr static double Q{0.0};
    constexpr static double k{0.1};
    constexpr static double h = 0.1;
    constexpr static double U0 = 1.0;
    static double alphaImpl() {
        double Pe = U0*h/(2*k);
        if(Pe == 0) {
            return 0.0;
        }
        else {
            return std::cosh(Pe)/std::sinh(Pe) - 1.0/Pe;
        }
    }

    // Pe = Uh/(2*k) < 1

    static Tensor<NSD, 1, double> advectionVelocity(coordinate<> const &x)
    {
        Tensor<NSD, 1, double> ret(0);
        ret(0) = U0;
        return ret;
    };
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_SUPG_ADVECTION_HPP
