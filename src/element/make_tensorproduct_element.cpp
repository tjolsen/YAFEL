//
// Created by tyler on 4/26/17.
//

#include "element/Element.hpp"
#include "element/ShapeFunctionUtils.hpp"

YAFEL_NAMESPACE_OPEN


void Element::make_tensorProduct()
{

    //auto lob_pts_1d = make_LobattoPoints_1D(elementType.polyOrder);
    auto qr = QuadratureRule(elementType.polyOrder + 1,
                             QuadratureRule::QuadratureType::GAUSS_LOBATTO);
    std::vector<double> lob_pts_1d;
    for (auto x : qr.nodes)
        lob_pts_1d.push_back(x(0));

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
        cellTypes.resize((npts - 1) * (npts - 1) * (npts - 1), CellType::Hex8);

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
        int j_stride = npts;
        int i_stride = npts * npts;
        for (int i = 0; i < npts - 1; ++i) {
            for (int j = 0; j < npts - 1; ++j) {
                for (int k = 0; k < npts - 1; ++k) {
                    int corner = i * i_stride + j * j_stride + k;

                    cellNodes[offset + 0] = corner;
                    cellNodes[offset + 1] = corner + 1;
                    cellNodes[offset + 2] = corner + 1 + j_stride;
                    cellNodes[offset + 3] = corner + j_stride;
                    cellNodes[offset + 4] = corner + i_stride;
                    cellNodes[offset + 5] = corner + i_stride + 1;
                    cellNodes[offset + 6] = corner + i_stride + 1 + j_stride;
                    cellNodes[offset + 7] = corner + i_stride + j_stride;

                    offsets[idx++] = offset;
                    offset += 8;

                }
            }
        }
        offsets[idx] = offset;

    } else {
        throw "Unsupported topoDim";
        // unsupported topoDim
    }

    // Set local mesh members from created containers
    localMesh.setGeometryNodes(localPoints_xi);
    localMesh.setCellNodes(cellNodes);
    localMesh.setOffsets(offsets);
    localMesh.setCellTypes(cellTypes);

    quadratureRule = QuadratureRule::make_tensor_product(QuadratureRule::QuadratureType::GAUSS_LEGENDRE,
                                                         elementType.topoDim, 2 * elementType.polyOrder);

    tensor_product_shape_functions(localMesh.getGeometryNodes(),
                                   quadratureRule.nodes,
                                   elementType.topoDim,
                                   shapeValues,
                                   shapeGradXi);

    if(elementType.topoDim == 2) {
        build_quad_faces();
    }
}



YAFEL_NAMESPACE_CLOSE