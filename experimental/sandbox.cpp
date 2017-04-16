//
// Created by tyler on 3/14/17.
//

#include "mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include <element/Element.hpp>
#include <element/ElementFactory.hpp>
#include "utils/DoFManager.hpp"
#include "output/VTUBackend.hpp"
#include "output/OutputMesh.hpp"
#include "output/OutputFrame.hpp"
#include "boundary_conditions/DirichletBC.hpp"

#include <iostream>

using namespace yafel;

int main()
{
    Mesh M("minsquare.msh");
    int p = 3;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, 1);

    auto f = [](const coordinate<>&x,double t){return x(0) + t;};

    DirichletBC bc(dofm, [](const coordinate<>&x, double t){return x(0)+t;});

    bc.selectByFunction([](const coordinate<>&x) { return std::abs(x(0)) < 1.0e-8; }, 0);
    bc.selectByRegionID(1,0);
    bc.selectByRegionID(2,0);

    auto mask = bc.getMask();

    Eigen::VectorXd isBoundary(mask.size());
    for(auto n : IRange(0,static_cast<int>(mask.size()))) {
        isBoundary(n) = mask[n];
    }

    OutputMesh om(dofm);
    OutputFrame frame(om);
    OutputData data(isBoundary,"Is Boundary");
    frame.point_data.push_back(&data);

    VTUBackend vtu;
    vtu.initialize("boundaryOut",0);
    vtu.write_frame(frame);
    vtu.finalize();

    return 0;
}


