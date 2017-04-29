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

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>

#ifndef VIENNACL_WITH_EIGEN
#define VIENNACL_WITH_EIGEN 1
#endif

#include <viennacl/linalg/cg.hpp>
#include <iostream>



/**
 * \file
 *
 * Implementation of a basic steady-state linear elasticity solver
 * using the CG assembly framework. Includes use of a basic timer
 * (yafel::BasicTimer), which provides tic() and toc() functionality
 * to time the assembly process.
 *
 * The linear solver uses a ViennaCL OpenCL Conjugate Gradient solver,
 * offering drastic speed gains over a CPU-only (i.e. Eigen-only) solver
 * on my machine (which has an NVidia GTX 1070)
 */


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
        return;
        for (auto A : IRange(0, static_cast<int>(R_el.rows()))) {
            R_el(A) += 0 * E.shapeValues[qpi](A / NSD) * E.jxw;
        }
    }

    static void LocalTangent(const Element &E, int qpi, double,
                             Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &K_el)
    {

        //note permuted indices in tensor functor. this is so that the "j" and "l"
        // indices are adjacent in memory, which is necessary for the later operations
        Tensor<NSD, 4> Cikjl = make_TensorFunctor<NSD, 4, double>(
                [](int i, int k, int j, int l) {
                    auto lambda = LinearElasticity<NSD>::lambda;
                    auto mu = LinearElasticity<NSD>::mu;
                    return lambda * (i == j) * (k == l) + mu * ((i == k) * (j == l) + (i == l) * (j == k));
                }
        );


        const int nNodes = K_el.rows() / NSD;
        //loops attempting strength reduction by inserting more loops!
        int A{0};
        for (auto Anode : IRange(0, nNodes)) {
            auto grad_Aj = make_TensorMap<NSD, 1>(&E.shapeGrad(Anode, 0));

            for (auto i : IRange(0, NSD)) {
                int B{Anode * NSD};

                //Diagonal block
                auto tmp = otimes(grad_Aj, grad_Aj).eval();
                for (auto k : IRange(0, NSD)) {
                    K_el(A, B) -= dot(Cikjl(i, k, colon(), colon()), tmp) * E.jxw;
                    ++B;
                }

                //Off-diagonal blocks
                for (auto Bnode : IRange(Anode + 1, nNodes)) {
                    auto grad_Bl = make_TensorMap<NSD, 1>(&E.shapeGrad(Bnode, 0));
                    auto tmp = otimes(grad_Aj, grad_Bl).eval();
                    for (auto k : IRange(0, NSD)) {
                        K_el(A, B) -= dot(Cikjl(i, k, colon(), colon()), tmp) * E.jxw;
                        K_el(B, A) = K_el(A, B);

                        ++B;
                    }
                }

                ++A;
            }
        }

    }
};

Eigen::VectorXd solveSystem(const Eigen::SparseMatrix<double, Eigen::RowMajor> &A, const Eigen::VectorXd &rhs)
{
    std::cout << "Starting solve...\n" << std::endl;

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

    return result;


    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::Upper | Eigen::Lower> solver;
    //solver.compute(A);
    //return solver.solve(rhs);
}

int main()
{

    constexpr int NSD = 3;

    Mesh M("mesh.msh");
    int p = 2;
    int dofpn = NSD;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, dofpn);
    FESystem feSystem(dofm);
    BasicTimer timer;

    timer.tic();
    CGAssembly<LinearElasticity<NSD>>(feSystem);
    timer.toc();

    std::cout << "Assembly time: " << timer.duration<>() << " ms" << std::endl;

    DirichletBC bc00(dofm, 0.0, 0);
    bc00.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });
    DirichletBC bc01(dofm, 0.0, 1);
    bc01.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });
    DirichletBC bc02(dofm, 0.0, 2);
    bc01.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });


    DirichletBC bc10(dofm, 0.001, 0);
    bc10.selectByFunction([](const coordinate<> &x) { return std::abs(x(0) - 2) < 1.0e-6; });

    bc00.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc01.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc02.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc10.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    Eigen::VectorXd U(feSystem.getGlobalResidual().rows());
    U = solveSystem(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    OutputMesh outputMesh(dofm);
    OutputData data(
            U,
            "Solution",
            OutputData::DataLocation::Point,
            OutputData::DataType::Vector,
            std::vector<int>(NSD,1)
    );
    OutputFrame frame(outputMesh);
    frame.point_data.push_back(&data);

    VTUBackend vtuBackend;
    vtuBackend.initialize("basic_elasticity");
    vtuBackend.write_frame(frame);
    vtuBackend.finalize();

    return 0;
}