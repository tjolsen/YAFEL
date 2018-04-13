//
// Created by tyler on 4/10/17.
//

#include "mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include "element/Element.hpp"
#include "element/ElementFactory.hpp"
#include "utils/DoFManager.hpp"
#include <iostream>

using namespace yafel;

/**
 * Compute the area and perimeter of a 2D mesh given in 2x2square.msh (created by gmsh from 2x2square.geo)
 * After generating the mesh (in terminal: gmsh -2 2x2square.geo), this program should
 * produce the output "A = 4" and "Perimeter = 8" for any choice of polynomial interpolation order "p".
 * @return
 */

int main()
{
    constexpr int NSD = 3;

    Mesh M("2x2square.msh");
    M.buildInternalFaces();

    int p = 1; // set polynomial interpolation order ( p <= 5 for now if using triangle elements )
    DoFManager dofm(M, DoFManager::ManagerType::DG, p);
    double area{0};

    ElementFactory EF(1,1);

    for (auto c : IRange(0, M.nCells())) {
        auto et = dofm.CellType_to_ElementType(M.getCellType(c), p);

	if(et.topoDim != NSD) {
	    continue;
	}
	
        auto &E = EF.getElement(et);
        for (auto qpi : IRange(0, E.nQP())) {
            E.update<NSD>(c, qpi, dofm);
            area += E.jxw;
        }
    }
    std::cout << "A = " << area << std::endl;

    // Surface area (i.e. perimeter) calculation
    double Perimeter{0};
    std::vector<int> local_fnodes;
    for (int fi = 0; fi < M.getInternalFaces().size(); ++fi) {
        auto F = M.getInternalFaces()[fi];
        if (!(F.left < 0 || F.right < 0)) {
            //internal face
            continue;
        }
        int elem{-1};
        if (F.left >= 0) {
            elem = F.left;
            dofm.getLeftFaceNodes(fi, local_fnodes);
        }

        if (F.right >= 0) {
            elem = F.right;
            dofm.getRightFaceNodes(fi,local_fnodes);
        }
        auto et = dofm.element_types[elem];
        if(et.topoDim != NSD) {
            //this "internal face" is likely actually an internal edge.
            continue;
        }

        auto &E = EF.getElement(et);

        for (auto fqpi : IRange(0, E.nFQP())) {
            E.face_update<NSD>(elem, fqpi, local_fnodes, dofm);

            Perimeter += E.jxw;
        }
    }

    std::cout << "Perimeter = " << Perimeter << std::endl;

}
