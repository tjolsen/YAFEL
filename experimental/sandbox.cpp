//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "assembly/CGAssembly.hpp"
#include "boundary_conditions/DirichletBC.hpp"
#include "output/SimulationOutput.hpp"

#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>


using namespace yafel;

int main()
{
    Mesh M("minsquare.msh");
    M.buildInternalFaces();

    int p = 3;

    auto f = [](coordinate<> x, double t) {
        using std::sin;
        using std::cos;

        return cos(t)*sin(10*x(0))*cos(10*x(1));
    };

    DoFManager dofm(M, DoFManager::ManagerType::DG, p);
    FESystem feSystem(dofm);
    auto &U = feSystem.getSolution();
    std::vector<int> container;
    SimulationOutput simulationOutput("output", BackendType::HDF5);
    double dt = 0.1;

    for(int ti=0; ti < 100; ++ti) {
        double t = dt*ti;
        feSystem.currentTime() = t;

        for (auto e : IRange(0, dofm.nCells())) {
            dofm.getGlobalNodes(e, container);
            for (auto n : container) {
                auto x = dofm.dof_nodes[n];
                U(n) = f(x,t);
            }
        }

        simulationOutput.captureFrame(feSystem);
        int a = 1;
    }



    return 0;
}


