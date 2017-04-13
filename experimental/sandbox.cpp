//
// Created by tyler on 3/14/17.
//

#include "mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include <iostream>
#include <element/Element.hpp>
#include <element/ElementFactory.hpp>
#include "utils/DoFManager.hpp"

#include "output/VTUBackend.hpp"
#include "output/OutputMesh.hpp"
#include "output/OutputFrame.hpp"

using namespace yafel;

int main()
{
    Mesh M("singleTet.msh");
    int p = 4;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, 1);
    ElementFactory EF(1);
    auto &E = EF.getElement(dofm.element_types[0]);


    for(auto x : E.localMesh.getGeometryNodes()) {
        std::cout << x(0) << "  " << x(1) << "  " << x(2) << std::endl;
    }

    return 0;
}


