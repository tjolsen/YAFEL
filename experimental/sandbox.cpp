//
// Created by tyler on 3/14/17.
//

#include "new_mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include <iostream>
#include <element/Element.hpp>
#include <element/ElementFactory.hpp>
#include "utils/DoFManager.hpp"

using namespace yafel;

int main()
{
    /*
    //---------------------------------------------------
    std::vector<coordinate<>> X;
    int N = 3;
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
    */
    Mesh M("minsquare.msh");
    int p=1;
    DoFManager dofm(M,DoFManager::ManagerType::CG, p, 1);

    double area{0};
    ElementFactory EF(1);

    for(auto c : IRange(0,M.nCells())) {
        auto et = dofm.CellType_to_ElementType(M.getCellType(c),p);
        auto &E = EF.getElement(et);

        for(auto qpi : IRange(0,E.nQP())) {
            E.update(qpi, dofm);

            auto w = E.quadratureRule.weights[qpi];
            area += w*E.detJ;
        }
    }


    return 0;
}


