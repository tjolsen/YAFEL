//
// Created by tyler on 4/17/17.
//

#include "yafel_globals.hpp"
#include "assembly/CGAssembly.hpp"
#include "boundary_conditions/DirichletBC.hpp"
#include "output/SimulationOutput.hpp"
#include "utils/BasicTimer.hpp"

#ifndef VIENNACL_WITH_EIGEN
#define VIENNACL_WITH_EIGEN 1
#endif

#include <viennacl/linalg/cg.hpp>
#include <iostream>
#include <eigen3/Eigen/IterativeLinearSolvers>

using namespace yafel;

template<int NSD>
struct PoissonEquation
{
    static constexpr int nsd() { return NSD; }

    static void
    LocalResidual(const Element &E, int qpi, double,
                  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &R_el)
    {

        R_el += (100.0*E.jxw) * E.shapeValues[qpi];
        /*
        for (auto A : IRange(0, static_cast<int>(R_el.rows()))) {
            R_el(A) = 10*E.jxw*E.shapeValues[qpi](A);
        }*/

    }

    static void LocalTangent(const Element &E, int, double,
                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &K_el)
    {
        K_el -= E.shapeGrad * (E.shapeGrad.transpose() * E.jxw);
        /*
        for (auto A: IRange(0, static_cast<int>(K_el.rows()))) {
            auto grad_Ai = make_TensorMap<NSD,1>(&E.shapeGrad(A,0));
            for (auto B: IRange(A, static_cast<int>(K_el.rows()))) {
                auto grad_Bi = make_TensorMap<NSD,1>(&E.shapeGrad(B,0));
                K_el(A, B) -=  dot(grad_Ai, grad_Bi)* E.jxw;
                K_el(B, A) = K_el(A, B);
            }

        }*/
    }
};

Eigen::VectorXd solveSystem(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &rhs)
{

    std::cout << "Norm0 = " << rhs.norm() << std::endl;


    viennacl::vector<double> vcl_rhs(rhs.rows());
    viennacl::compressed_matrix<double> vcl_A(A.rows(), A.cols());
    viennacl::copy(rhs, vcl_rhs);
    viennacl::copy(A, vcl_A);

    //solve
    viennacl::linalg::cg_tag tag(1.0e-14, 2 * rhs.rows());
    viennacl::vector<double> vcl_result = viennacl::linalg::solve(vcl_A, vcl_rhs, tag);

    std::cout << ": System solved.\n\t Iterations: " << tag.iters() << "\n\tError: " << tag.error() << std::endl;

    //Copy back
    Eigen::VectorXd result(rhs.rows());
    viennacl::copy(vcl_result, result);


    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper | Eigen::Lower> solver;
    //solver.compute(A);
    //Eigen::VectorXd result = solver.solve(rhs);
    std::cout << "new norm = " << (A * result - rhs).norm() << std::endl;

    return result;
}

int main()
{
    Mesh M("mesh.msh");
    int p = 2;
    int dofpn = 1;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, dofpn);
    FESystem feSystem(dofm);
    BasicTimer timer;
    timer.tic();
    CGAssembly<PoissonEquation<3>>(feSystem);
    timer.toc();

    std::cout << "Assembly time: " << timer.duration<std::chrono::microseconds>() << " us" << std::endl;

    DirichletBC bc0(dofm, 0.0);
    bc0.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });

    DirichletBC bc1(dofm, 0.0);
    bc1.selectByFunction([](const coordinate<> &x) { return std::abs(x(0) - 2) < 1.0e-6; });

    bc0.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc1.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    auto &U = feSystem.getSolution();
    U = solveSystem(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());


    SimulationOutput simulationOutput("output", BackendType::HDF5);
    simulationOutput.captureFrame(feSystem);

    return 0;
}