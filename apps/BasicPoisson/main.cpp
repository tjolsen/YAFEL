//
// Created by tyler on 4/17/17.
//

#include "yafel_globals.hpp"
#include "assembly/CGAssembly.hpp"
#include "boundary_conditions/DirichletBC.hpp"

#include "output/OutputData.hpp"
#include "output/OutputMesh.hpp"
#include "output/OutputFrame.hpp"
#include "output/VTUBackend.hpp"

#include "utils/BasicTimer.hpp"

#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>


using namespace yafel;

template<int NSD>
struct PoissonEquation
{
    static constexpr int nsd() { return NSD; }

    static void
    LocalResidual(const Element &E, int qpi, double, Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &R_el)
    {
        for (auto A : IRange(0, static_cast<int>(R_el.rows()))) {
            R_el(A) += 10*E.shapeValues[qpi](A) * E.jxw;
        }
    }

    static void LocalTangent(const Element &E, int, double,
                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &K_el)
    {
        for (auto A: IRange(0, static_cast<int>(K_el.rows()))) {
            auto grad_Ai = make_TensorMap<NSD,1>(&E.shapeGrad(A,0));
            for (auto B: IRange(A, static_cast<int>(K_el.rows()))) {
                auto grad_Bi = make_TensorMap<NSD,1>(&E.shapeGrad(A,0));
                K_el(A, B) -=  dot(grad_Ai, grad_Bi)* E.jxw;
                K_el(B, A) = K_el(A, B);
            }

        }
    }
};

Eigen::VectorXd solveSystem(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &rhs)
{
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper | Eigen::Lower> solver;
    solver.compute(A);
    return solver.solve(rhs);
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
    CGAssembly<PoissonEquation<2>>(feSystem);
    timer.toc();

    std::cout << "Assembly time: " << timer.duration<std::chrono::microseconds>() << " us" << std::endl;

    DirichletBC bc0(dofm, [](const coordinate<> &x,double){return x(0)*x(0)/4;});
    bc0.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6 || std::abs(x(1)) < 1.0e-6; });

    DirichletBC bc1(dofm, 1.0);
    bc1.selectByFunction([](const coordinate<> &x) { return std::abs(x(0) - 2) < 1.0e-6; });

    bc0.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc1.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    Eigen::VectorXd U(feSystem.getGlobalResidual().rows());
    U = solveSystem(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    OutputMesh outputMesh(dofm);
    OutputData data(U, "Solution");
    OutputFrame frame(outputMesh);
    frame.point_data.push_back(&data);

    VTUBackend vtuBackend;
    vtuBackend.initialize("basic_poisson");
    vtuBackend.write_frame(frame);
    vtuBackend.finalize();

    return 0;
}