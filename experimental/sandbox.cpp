//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "assembly/CGAssembly.hpp"
#include "element/element_boundary_nodes.hpp"
#include "boundary_conditions/DirichletBC.hpp"

#include "output/SimulationOutput.hpp"
//#include "output/OutputData.hpp"
//#include "output/OutputMesh.hpp"
//#include "output/OutputFrame.hpp"
//#include "output/VTUBackend.hpp"

#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>



using namespace yafel;

int main()
{
    Mesh M("minsquare.msh");
    M.buildInternalFaces();

    int p = 4;

    DoFManager dofm(M, DoFManager::ManagerType::DG, p);
    FESystem feSystem(dofm);
    auto &U = feSystem.getSolution();
    std::vector<int> container;
    for(auto e : IRange(0,dofm.nCells())) {
        dofm.getGlobalNodes(e,container);
        for(auto n : container) {
            U(n) = e;
        }
    }

    SimulationOutput simulationOutput("output", BackendType::HDF5);

    simulationOutput.captureFrame(feSystem);

    return 0;
}


