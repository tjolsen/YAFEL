//
// Created by tyler on 4/17/17.
//

#include "yafel_globals.hpp"
#include "assembly/CGAssembly.hpp"
#include "assembly/ZZGradientRecovery.hpp"
#include "boundary_conditions/DirichletBC.hpp"
#include "output/SimulationOutput.hpp"
#include "utils/BasicTimer.hpp"
#include "lin_alg/linear_solvers/LinearSolve.hpp"

/*
#ifndef VIENNACL_WITH_EIGEN
#define VIENNACL_WITH_EIGEN 1
#endif

#include <viennacl/linalg/cg.hpp>
#include <iostream>
#include <Eigen/IterativeLinearSolvers>
*/

using namespace yafel;

template<int NSD>
struct PoissonEquation
{
    static constexpr int nsd()
    { return NSD; }

    static void
    LocalResidual(const Element &E, int qpi, coordinate<>&, double,
                  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &u_el,
                  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &R_el)
    {

        R_el += (100.0 * E.jxw) * E.shapeValues[qpi];
        /*
        for (auto A : IRange(0, static_cast<int>(R_el.rows()))) {
            R_el(A) = 10*E.jxw*E.shapeValues[qpi](A);
        }*/

    }

    static void LocalTangent(const Element &E, int qpi, coordinate<>&, double,
                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &u_el,
                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &K_el)
    {
        //K_el -= E.shapeGrad* (E.shapeGrad.transpose() * E.jxw);

        for (auto A: IRange(0, static_cast<int>(K_el.rows()))) {
            auto grad_Ai = make_TensorMap<NSD,1>(&E.shapeGrad(A,0));
            for (auto B: IRange(A, static_cast<int>(K_el.rows()))) {
                auto grad_Bi = make_TensorMap<NSD,1>(&E.shapeGrad(B,0));
                K_el(A, B) -=  dot(grad_Ai, grad_Bi)* E.jxw;
                K_el(B, A) = K_el(A, B);
            }

        }
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


    /*
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper | Eigen::Lower> solver;
    solver.compute(A);
    Eigen::VectorXd result = solver.solve(rhs);
    std::cout << "new norm = " << (A * result - rhs).norm() << std::endl;
     */

    return result;
}

int main()
{
    constexpr int nsd = 3;
    Mesh M("mesh.msh");
    int p = 1;
    int dofpn = 1;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, dofpn);
    FESystem feSystem(dofm);
    BasicTimer timer;
    timer.tic();
    CGAssembly<PoissonEquation<nsd>>(feSystem);
    timer.toc();

    std::cout << "Assembly time: " << timer.duration<std::chrono::microseconds>() << " us" << std::endl;

    DirichletBC bc0(dofm, 0.0);
    bc0.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });

    //DirichletBC bc1(dofm, 0.0);
    //bc1.selectByFunction([](const coordinate<> &x) { return std::abs(x(1)) < 1.0e-6; });

    bc0.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    //bc1.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    auto &U = feSystem.getSolution();
    U = solveSystem(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    std::cout << "Nnodes = " << dofm.nNodes() << std::endl;

    ZZGradientRecovery<PoissonEquation<nsd>>(feSystem);

    std::function<void(FESystem&,OutputFrame&)> captureFunc = [](FESystem &feSys, OutputFrame &frame) {
        frame.time = 0;
        OutputData::DataLocation dataLocation = OutputData::DataLocation::Point;
        std::string sol_name = "U";
        OutputData::DataType dt = OutputData::DataType::Scalar;
        auto dat = std::make_shared<OutputData>(feSys.getSolution(),
                                                sol_name, dataLocation, dt, std::vector<int>{1});


        frame.addData(dat);

        int nsd = feSys.getSolutionGradient().cols();
        auto & gradient = feSys.getSolutionGradient();

        auto grad_dat = std::make_shared<OutputData>(
                gradient,
                "gradU",
                OutputData::DataLocation::Point,
                OutputData::DataType::Vector,
                std::vector<int>(nsd, 1));

        frame.addData(grad_dat);
    };

    SimulationOutput simulationOutput("output", BackendType::VTU);
    simulationOutput.captureFrame(feSystem, captureFunc);

    //std::cout << feSystem.getGlobalTangent() << std::endl << std::endl << feSystem.getGlobalResidual() << std::endl;

    return 0;
}
