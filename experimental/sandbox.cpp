//
// Created by tyler on 3/14/17.
//

#include "element/Element.hpp"
#include <iostream>

using namespace yafel;


int main() {

    // make a quad
    Element e(ElementType(ElementClass::TensorProduct,3,3));


    for(auto x : e.localMesh.getGeometryNodes())
    {
        std::cout << x(0) << "  " << x(1) << "  " << x(2) << std::endl;
    }
    std::vector<int> cellNodes(2);
    for(int c=0; c<e.localMesh.nCells(); ++c) {
        e.localMesh.getCellNodes(c,cellNodes);
        std::cout << c << ":  ";
        for(auto n : cellNodes) {
            std::cout << n << "   ";
        }
        std::cout << std::endl;
    }

    return 0;
}


