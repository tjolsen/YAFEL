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

    static void LocalTangent(const Element &E, int, double,
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
        for(int Anode=0; Anode < nNodes; ++Anode) {
            auto grad_Aj = make_TensorMap<NSD, 1>(&E.shapeGrad(Anode, 0));

            for(int i=0; i<NSD; ++i) {
                int B{Anode * NSD};

                //Diagonal block
                auto tmp = otimes(grad_Aj, grad_Aj).eval();
                for(int k=0; k<NSD; ++k){
                    K_el(A, B) -= dot(Cikjl(i, k, colon(), colon()).eval(), tmp) * E.jxw;
                    ++B;
                }

                //Off-diagonal blocks
                for(int Bnode = Anode+1; Bnode<nNodes; ++Bnode) {
                    auto grad_Bl = make_TensorMap<NSD, 1>(&E.shapeGrad(Bnode, 0));
                    auto tmp = otimes(grad_Aj, grad_Bl).eval();
                    //for (auto k : IRange(0, NSD)) {
                    for(int k=0; k<NSD; ++k) {
                        K_el(A, B) -= dot(Cikjl(i, k, colon(), colon()).eval(), tmp) * E.jxw;
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
    BasicTimer timer;

    timer.tic();
    Mesh M("mesh.msh");
    timer.toc();

    std::cout <<"Mesh import time: " << timer.duration<>() << " ms" << std::endl;
    int p = 1;
    int dofpn = NSD;

    timer.tic();
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, dofpn);
    timer.toc();
    std::cout <<"dofm creation time: " << timer.duration<>() << " ms" << std::endl;

    FESystem feSystem(dofm);


    timer.tic();
    CGAssembly<LinearElasticity<NSD>>(feSystem);
    timer.toc();

    std::cout << "Assembly time: " << timer.duration<>() << " ms" << std::endl;


    DirichletBC bc00(dofm, 0.0, 0);
    bc00.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });
    DirichletBC bc01(dofm, 0.0, 1);
    bc01.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });
    DirichletBC bc02(dofm, 0.0, 2);
    bc02.selectByFunction([](const coordinate<> &x) { return std::abs(x(0)) < 1.0e-6; });

    double twist_angle = 0.1; // angle in radians
    auto twist = [twist_angle](auto x) {
        x(0) = 0;
        Tensor<3,1> dx = x - Tensor<3,1>{0,1,1};
        auto r = norm(dx);
        auto theta = std::atan2(dx(2), dx(1));

        return (twist_angle*Tensor<3,1>(0, -r*std::sin(theta), r*std::cos(theta))).eval();
    };

    DirichletBC bc10(dofm, 0.0, 0);
    DirichletBC bc11(dofm, [twist](auto x, double){return twist(x)(1);}, 1);
    DirichletBC bc12(dofm, [twist](auto x, double){return twist(x)(2);}, 2);
    bc10.selectByFunction([](const coordinate<> &x) { return std::abs(x(0) - 2) < 1.0e-6; });
    bc11.selectByFunction([](const coordinate<> &x) { return std::abs(x(0) - 2) < 1.0e-6; });
    bc12.selectByFunction([](const coordinate<> &x) { return std::abs(x(0) - 2) < 1.0e-6; });

    bc00.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc01.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc02.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc10.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc11.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc12.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    timer.tic();
    Eigen::VectorXd U(feSystem.getGlobalResidual().rows());
    U = solveSystem(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    timer.toc();
    std::cout <<"Solution time: " << timer.duration<>() << " ms" << std::endl;

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