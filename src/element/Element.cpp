//
// Created by tyler on 3/14/17.
//

#include <new_mesh/LobattoPoints1D.hpp>
#include "element/Element.hpp"

YAFEL_NAMESPACE_OPEN

Element::Element(ElementType et)
        : elementType(et), localMesh(Mesh::DefinitionScheme::Explicit)
{

    switch (et.elementClass) {
        case ElementClass::TensorProduct:
            make_tensorProduct();
            break;
        case ElementClass::Simplex:
            make_simplex();
            break;
        case ElementClass::None:
            //null element likely used for
            break;
    }
}


void Element::make_tensorProduct()
{

    auto lob_pts_1d = make_LobattoPoints_1D(elementType.polyOrder);
    auto npts = static_cast<int>(lob_pts_1d.size());

    std::vector<coordinate<>> localPoints_xi;
    std::vector<int> cellNodes;
    std::vector<int> offsets;
    std::vector<CellType> cellTypes;

    if (elementType.topoDim == 1) {
        //make a line element
        localPoints_xi.resize(npts);
        cellNodes.resize(2 * (npts - 1));
        offsets.resize(npts);
        cellTypes.resize(npts - 1, CellType::Line2);


        for (int i = 0; i < npts; ++i) {
            localPoints_xi[i](0) = lob_pts_1d[i];
        }

        int offset = 0;
        for (int i = 0; i < npts - 1; ++i) {
            offsets[i] = offset;
            offset += 2;

            cellNodes[2 * i] = i;
            cellNodes[2 * i + 1] = i + 1;
        }
        offsets[npts - 1] = offset;

    } else if (elementType.topoDim == 2) {
        // make a quad element

        localPoints_xi.resize(npts * npts);
        cellNodes.resize(4 * (npts - 1) * (npts - 1));
        offsets.resize((npts - 1) * (npts - 1) + 1);
        cellTypes.resize((npts - 1) * (npts - 1), CellType::Quad4);

        int idx = 0;
        for (int i = 0; i < npts; ++i) {
            for (int j = 0; j < npts; ++j) {
                localPoints_xi[idx](1) = lob_pts_1d[i];
                localPoints_xi[idx](0) = lob_pts_1d[j];
                ++idx;
            }
        }

        int offset = 0;
        idx = 0;
        for (int i = 0; i < npts - 1; ++i) {
            for (int j = 0; j < npts - 1; ++j) {

                cellNodes[offset + 0] = i * npts + j;
                cellNodes[offset + 1] = i * npts + j + 1;
                cellNodes[offset + 2] = (i + 1) * npts + j + 1;
                cellNodes[offset + 3] = (i + 1) * npts + j;

                offsets[idx++] = offset;
                offset += 4;
            }
        }
        offsets[idx] = offset;


    } else if (elementType.topoDim == 3) {
        // make a hex element

        localPoints_xi.resize(npts * npts * npts);
        cellNodes.resize(8 * (npts - 1) * (npts - 1) * (npts - 1));
        offsets.resize((npts - 1) * (npts - 1) * (npts - 1) + 1);
        cellTypes.resize((npts - 1) * (npts - 1) * (npts - 1), CellType::Quad4);

        int idx = 0;
        for (int i = 0; i < npts; ++i) {
            for (int j = 0; j < npts; ++j) {
                for (int k = 0; k < npts; ++k) {
                    localPoints_xi[idx](2) = lob_pts_1d[i];
                    localPoints_xi[idx](1) = lob_pts_1d[j];
                    localPoints_xi[idx](0) = lob_pts_1d[k];
                    ++idx;
                }
            }
        }

        int offset = 0;
        idx = 0;
        for (int i = 0; i < npts - 1; ++i) {
            for (int j = 0; j < npts - 1; ++j) {
                for (int k = 0; k < npts - 1; ++k) {

                    cellNodes[offset + 0] = (i * npts + j) * npts + k;
                    cellNodes[offset + 1] = (i * npts + j) * npts + k + 1;
                    cellNodes[offset + 2] = (i * npts + j + 1) * npts + k + 1;
                    cellNodes[offset + 3] = (i * npts + j + 1) * npts + k;
                    cellNodes[offset + 4] = ((i + 1) * npts + j) * npts + k;
                    cellNodes[offset + 5] = ((i + 1) * npts + j) * npts + k + 1;
                    cellNodes[offset + 6] = ((i + 1) * npts + j + 1) * npts + k + 1;
                    cellNodes[offset + 7] = ((i + 1) * npts + j + 1) * npts + k;

                    offsets[idx++] = offset;
                    offset += 8;


                }
            }
        }

    } else {
        // unsupported topoDim

    }



    // Set local mesh members from created containers
    localMesh.setGeometryNodes(std::move(localPoints_xi));
    localMesh.setCellNodes(std::move(cellNodes));
    localMesh.setOffsets(std::move(offsets));
    localMesh.setCellTypes(std::move(cellTypes));

}


void Element::make_simplex()
{}

YAFEL_NAMESPACE_CLOSE