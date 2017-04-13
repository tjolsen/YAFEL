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
    Mesh M("minsquare.msh");
    int p = 1;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, 1);

    M.buildInternalFaces();

    for(auto i : M.getInternalFaces()) {
        std::cout << i << std::endl;
    }


    std::cout << std::endl << std::endl;
    for(auto i : M.getBoundaryFaceIdxs()) {
        std::cout << i << "  "
                  << M.getInternalFaces()[i]
                  << std::endl;
    }




    return 0;
}


