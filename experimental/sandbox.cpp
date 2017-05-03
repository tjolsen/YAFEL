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

    int p = 2;

    auto f = [](coordinate<> x) { return x(1); };

    DoFManager dofm(M, DoFManager::ManagerType::DG, p);
    FESystem feSystem(dofm);
    auto &U = feSystem.getSolution();
    std::vector<int> container;
    SimulationOutput simulationOutput("output", BackendType::HDF5);

    for (auto e : IRange(0, dofm.nCells())) {
        dofm.getGlobalNodes(e, container);
        for (auto n : container) {
            auto x = dofm.dof_nodes[n];
            U(n) = f(x);
        }
    }

    ElementFactory EF(1);
    for(auto e : IRange(0,dofm.nCells())) {
        int qpi = 1;
        dofm.getGlobalNodes(e,container);
        auto &E = EF.getElement(dofm.element_types[e]);
        E.update<2>(e,qpi,dofm);

        Tensor<2,1> Ugrad(0);
        for(int i=0; i<container.size(); ++i) {

            auto gradmap = make_TensorMap<2,1>(&E.shapeGrad(i,0));

            Ugrad = Ugrad + U(container[i])*gradmap;
        }

        std::cout << "grad(U) = {" << Ugrad(0) << ", " << Ugrad(1) << "}\n";
    }


    simulationOutput.captureFrame(feSystem);



    return 0;
}


