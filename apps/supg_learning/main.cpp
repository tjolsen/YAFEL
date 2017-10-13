//
// Created by tyler on 10/12/17.
//

#include "yafel.hpp"
#include "supg_advection.hpp"

#include <eigen3/Eigen/IterativeLinearSolvers>

using namespace yafel;

template<typename MatType>
auto solveSystem(const MatType &A, const Eigen::VectorXd &b) {

    Eigen::BiCGSTAB<MatType> solver;
    solver.compute(A);
    Eigen::VectorXd x = solver.solve(b);
    return x;
}

int main() {


    constexpr int nsd = 2;
    Mesh M("mesh.msh");
    int p=2;
    int dofpn = 1;

    DoFManager dofm(M,DoFManager::ManagerType::CG, p, dofpn);
    FESystem feSystem(dofm);
    CGAssembly<AdvectionDiffusionSUPG<nsd>>(feSystem);

    DirichletBC bc0(dofm,0.0);
    bc0.selectByFunction([](auto x){return std::abs(x(0)) < 1.0e-6;});
    DirichletBC bc1(dofm,1.0);
    bc1.selectByFunction([](auto x){return std::abs(x(0) - 1) < 1.0e-6;});

    bc0.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
    bc1.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    //std::cout << feSystem.getGlobalTangent() << std::endl;
    feSystem.getSolution() = Eigen::VectorXd(feSystem.getGlobalResidual().rows());
    feSystem.getSolution() = solveSystem(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());

    SimulationOutput simulationOutput("output", BackendType::VTU);
    simulationOutput.captureFrame(feSystem);
}