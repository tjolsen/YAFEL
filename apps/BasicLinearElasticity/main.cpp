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
struct LinearElasticity
{
    static constexpr int nsd() { return NSD; }

    static constexpr double Youngs{200.0e9};
    static constexpr double nu{0.3};

    static constexpr double lambda = (Youngs * nu) / ((1 + nu) * (1 - 2 * nu));
    static constexpr double mu = Youngs / (2 * (1 + nu));

    static void
    LocalResidual(const Element &E, int qpi, double, Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &R_el)
    {
        for (auto A : IRange(0, static_cast<int>(R_el.rows()))) {
            R_el(A) += 0 * E.shapeValues[qpi](A / NSD) * E.jxw;
        }
    }

    static void LocalTangent(const Element &E, int qpi, double,
                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &K_el)
    {

        Tensor<NSD, 4> Cijkl = make_TensorFunctor<NSD, 4, double>(
                [](int i, int j, int k, int l) {
                    auto lambda = LinearElasticity<NSD>::lambda;
                    auto mu = LinearElasticity<NSD>::mu;
                    return lambda * (i == j) * (k == l) + mu * ((i == k) * (j == l) + (i == l) * (j == k));
                }
        );

        for (auto A: IRange(0, static_cast<int>(K_el.rows()))) {
            auto i = A % NSD;
            auto Anode = A / NSD;
            auto grad_Aj = make_TensorMap<NSD, 1>(&E.shapeGrad(Anode, 0));

            for (auto B: IRange(A, static_cast<int>(K_el.rows()))) {
                auto k = B % NSD;
                auto Bnode = B / NSD;
                auto grad_Bl = make_TensorMap<NSD, 1>(&E.shapeGrad(Bnode, 0));

                K_el(A, B) -= dot(grad_Aj, Cijkl(i, colon(), k, colon()) * grad_Bl);
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
    int p = 1;
    int dofpn = 2;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, dofpn);
    FESystem feSystem(dofm);
    CGAssembly<LinearElasticity<2>>(feSystem);


    DirichletBC bc00(dofm, 0.0, 0);
    bc00.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });
    DirichletBC bc01(dofm, 0.0, 1);
    bc01.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });

    DirichletBC bc10(dofm, 0.001, 0);
    bc10.selectByFunction([](const coordinate<> &x) { return std::abs(x(0) - 2) < 1.0e-6; });

    bc00.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc01.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc10.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    Eigen::VectorXd U(feSystem.getGlobalResidual().rows());
    U = solveSystem(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    OutputMesh outputMesh(dofm);
    OutputData data(
            U,
            "Solution",
            OutputData::DataLocation::Point,
            OutputData::DataType::Vector,
            {1, 1}
    );
    OutputFrame frame(outputMesh);
    frame.point_data.push_back(&data);

    VTUBackend vtuBackend;
    vtuBackend.initialize("basic_elasticity");
    vtuBackend.write_frame(frame);
    vtuBackend.finalize();

    return 0;
}