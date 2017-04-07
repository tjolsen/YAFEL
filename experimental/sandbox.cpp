//
// Created by tyler on 3/14/17.
//

#include "new_mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include <iostream>
#include <element/Element.hpp>

using namespace yafel;

int main()
{
    //---------------------------------------------------
    std::vector<coordinate<>> X;
    int N = 300;
    double dx = 1.0;
    for(auto i : IRange(0,N)) {
        for(auto j : IRange(0,N)) {
            coordinate<> x;
            x(0) = j*dx;
            x(1) = i*dx;
            X.push_back(x);
        }
    }

    std::vector<int> cells, cell_offsets;
    std::vector<CellType> celltypes;
    int offset = 0;
    for(auto i : IRange(0,N-1)) {
        for(auto j : IRange(0,N-1)) {
            int idx = i*N + j;
            cells.push_back(idx);
            cells.push_back(idx+1);
            cells.push_back(idx+1+N);
            cells.push_back(idx+N);
            cell_offsets.push_back(offset);
            offset += 4;
            celltypes.push_back(CellType::Quad4);
        }
    }
    cell_offsets.push_back(offset);

    //---------------------------------------------------

    Mesh M(Mesh::DefinitionScheme::Explicit,
           X,cells,cell_offsets,celltypes);


    M.buildInternalFaces();
    std::cout << M.getInternalFaces().size() << std::endl;
    auto &faces = M.getInternalFaces();
    for(int i=0; i<faces.size(); ++i) {
        std::cout << faces[i] << std::endl;
    }

    return 0;
}


