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
            R_el(A) += 0 * E.shapeValues[qpi](A) * E.jxw;
        }
    }

    static void LocalTangent(const Element &E, int qpi, double,
                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &K_el)
    {
        for (auto A: IRange(0, static_cast<int>(K_el.rows()))) {
            for (auto B: IRange(A, static_cast<int>(K_el.rows()))) {
                for (auto i : IRange(0, NSD)) {
                    K_el(A, B) -= E.shapeGrad(A, i) * E.shapeGrad(B, i) * E.jxw;
                }

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
    int p = 5;
    int dofpn = 1;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, dofpn);
    FESystem feSystem(dofm);
    CGAssembly<PoissonEquation<2>>(feSystem);

    DirichletBC bc0(dofm, 0.0);
    bc0.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });

    DirichletBC bc1(dofm, 1.0);
    bc1.selectByFunction([](const coordinate<> & x) { return std::abs(x(0) - 2) < 1.0e-6; });

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