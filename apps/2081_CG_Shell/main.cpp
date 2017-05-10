//
// Created by tyler on 5/9/17.
//

#include "yafel_globals.hpp"
#include "assembly/CGAssembly.hpp"

#include "output/OutputData.hpp"
#include "output/OutputMesh.hpp"
#include "output/OutputFrame.hpp"
#include "output/VTUBackend.hpp"
#include "utils/BasicTimer.hpp"

#include "LinearElasticity.hpp"
#include "SimplySupportedSquare.hpp"


#include <iostream>
#include <output/SimulationOutput.hpp>

using namespace yafel;

Eigen::VectorXd solveSystem(const Eigen::SparseMatrix<double, Eigen::RowMajor> &A, const Eigen::VectorXd &rhs)
{
    std::cout << "Starting solve...\n" << std::endl;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::Upper | Eigen::Lower> solver;
    solver.compute(A);
    std::cout << "Really Starting Solve..." << std::endl;
    return solver.solve(rhs);

}


int main(int argc, char **argv)
{

    constexpr int NSD = 3;
    int dofpn = NSD;
    std::string mesh_fname("thinPlate.msh");
    std::string output_fname("output");
    int polyOrder = 3;
    if (argc >= 4) {
        polyOrder = atoi(argv[1]);
        mesh_fname = argv[2];
        output_fname = argv[3];
    }

    output_fname += "_p_" + std::to_string(polyOrder);

    Mesh M(mesh_fname);
    DoFManager dofm(M, DoFManager::ManagerType::CG, polyOrder, dofpn);
    FESystem feSystem(dofm);

    CGAssembly<LinearElasticity<NSD>>(feSystem);
    {
        auto bcs = SimplySupportedSquare(dofm);
        for (auto &bc : bcs) {
            bc.apply(feSystem.getGlobalTangent(), feSystem.getGlobalResidual());
        }
    }

    feSystem.getSolution() = solveSystem(feSystem.getGlobalTangent(),feSystem.getGlobalResidual());

    SimulationOutput simulationOutput(output_fname, BackendType::HDF5);
    simulationOutput.captureFrame(feSystem);
    return 0;
}