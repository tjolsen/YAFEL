//
// Created by tyler on 3/14/17.
//

#include "element/Element.hpp"
#include "element/ShapeFunctionUtils.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

Element::Element(ElementType et, int dofpn)
        : elementType(et),
          localMesh(Mesh::DefinitionScheme::Explicit),
          dof_per_node(dofpn)
{

    switch (et.elementTopology) {
        case ElementTopology::TensorProduct:
            make_tensorProduct();
            break;
        case ElementTopology::Simplex:
            make_simplex();
            break;
        case ElementTopology::None:
            //null element. used for points/unsupported types
            break;
    }
}


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


void Element::make_simplex()
{
    auto qr = QuadratureRule(elementType.polyOrder + 1,
                             QuadratureRule::QuadratureType::GAUSS_LOBATTO);
    double snap_value = 1.0/(2<<20);
    auto snap = [snap_value](auto x)->coordinate<> { return round(x / snap_value) * snap_value; };

    std::vector<double> lob_pts_1d;
    for (auto x : qr.nodes)
        lob_pts_1d.push_back(x(0));
    auto npts = static_cast<int>(lob_pts_1d.size());

    std::vector<coordinate<>> localPoints_xi;
    std::vector<int> cellNodes;
    std::vector<int> offsets;
    std::vector<CellType> cellTypes;


    if (elementType.topoDim == 1) {
        //1D simplex is a line, which is also a 1D tensor-product element
        make_tensorProduct();
        return;
    } else if (elementType.topoDim == 2) {
        // based on scheme be Blyth & Pozrikdis (JAM 2005)
        localPoints_xi.resize((npts * (npts + 1)) / 2);
        int nt_local = elementType.polyOrder * elementType.polyOrder;
        cellNodes.resize(3 * nt_local);
        offsets.reserve(nt_local + 1);
        cellTypes.resize(nt_local, CellType::Tri3);

        int idx = 0;
        for (int i = 0; i < npts; ++i) {
            for (int j = 0; j < npts - i; ++j) {
                int k = npts - i - j - 1;
                double vi = (lob_pts_1d[i] + 1) / 2;
                double vj = (lob_pts_1d[j] + 1) / 2;
                double vk = (lob_pts_1d[k] + 1) / 2;
                localPoints_xi[idx](0) = (1 + 2 * vj - vi - vk) / 3;
                localPoints_xi[idx](1) = (1 + 2 * vi - vj - vk) / 3;
                ++idx;
            }
        }

        //make local triangles
        idx = 0;
        int offset = 0;
        for (int layer = 0; layer < elementType.polyOrder; ++layer) {
            for (int i = 0; i < elementType.polyOrder - layer; ++i) {

                cellNodes[offset + 0] = idx + i;
                cellNodes[offset + 1] = idx + i + 1;
                cellNodes[offset + 2] = idx + i + 1 + elementType.polyOrder - layer;
                offsets.push_back(offset);
                offset += 3;
                if (i < elementType.polyOrder - layer - 1) {
                    cellNodes[offset + 0] = idx + i + 1;
                    cellNodes[offset + 1] = idx + i + 2 + elementType.polyOrder - layer;
                    cellNodes[offset + 2] = idx + i + 1 + elementType.polyOrder - layer;
                    offsets.push_back(offset);
                    offset += 3;
                }

            }
            idx += elementType.polyOrder + 1 - layer;
        }
        offsets.push_back(offset);


    } else if (elementType.topoDim == 3) {
        // based on scheme in Blyth & Pozrikdis (JAM 2005)
        localPoints_xi.resize((npts * (npts + 1) * (npts + 2)) / 6);

        int idx{0};
        //xi-eta plane
        for (auto i : IRange(1, npts + 1)) {
            for (auto j : IRange(1, npts + 2 - i)) {
                int k = npts + 2 - i - j;
                double vi = (lob_pts_1d[i - 1] + 1) / 2;
                double vj = (lob_pts_1d[j - 1] + 1) / 2;
                double vk = (lob_pts_1d[k - 1] + 1) / 2;
                localPoints_xi[idx](0) = std::max(0.0, (1 + 2 * vj - vi - vk) / 3);
                localPoints_xi[idx](1) = std::max(0.0, (1 + 2 * vi - vj - vk) / 3);
                localPoints_xi[idx](2) = 0;
                ++idx;
            }
        }

        //eta-zeta plane
        for (auto j : IRange(1, npts)) {
            for (auto k : IRange(2, npts + 2 - j)) {
                int l = npts + 2 - j - k;
                double vj = (lob_pts_1d[j - 1] + 1) / 2;
                double vk = (lob_pts_1d[k - 1] + 1) / 2;
                double vl = (lob_pts_1d[l - 1] + 1) / 2;
                localPoints_xi[idx](0) = 0;
                localPoints_xi[idx](1) = std::max(0.0, (1 + 2 * vj - vk - vl) / 3);
                localPoints_xi[idx](2) = std::max(0.0, (1 + 2 * vk - vj - vl) / 3);
                ++idx;
            }
        }

        //xi-zeta plane
        for (auto i : IRange(2, npts)) {
            for (auto k : IRange(2, npts + 2 - i)) {
                int l = npts + 2 - i - k;
                double vi = (lob_pts_1d[i - 1] + 1) / 2;
                double vk = (lob_pts_1d[k - 1] + 1) / 2;
                double vl = (lob_pts_1d[l - 1] + 1) / 2;

                localPoints_xi[idx](0) = std::max(0.0, (1 + 2 * vi - vk - vl) / 3);
                localPoints_xi[idx](1) = 0;
                localPoints_xi[idx](2) = std::max(0.0, (1 + 2 * vk - vi - vl) / 3);

                ++idx;
            }
        }

        //slanted face
        for (auto i : IRange(2, npts)) {
            for (auto j : IRange(2, npts + 1 - i)) {
                int l = npts + 2 - i - j;
                double vi = (lob_pts_1d[i - 1] + 1) / 2;
                double vj = (lob_pts_1d[j - 1] + 1) / 2;
                double vl = (lob_pts_1d[l - 1] + 1) / 2;

                localPoints_xi[idx](0) = std::max(0.0, (1 + 2 * vi - vj - vl) / 3);
                localPoints_xi[idx](1) = std::max(0.0, (1 + 2 * vj - vi - vl) / 3);
                localPoints_xi[idx](2) = 1 - localPoints_xi[idx](0) - localPoints_xi[idx](1);

                ++idx;
            }
        }


        //interior points
        for (auto i : IRange(2, npts)) {
            for (auto j : IRange(2, npts + 1 - i)) {
                for (auto k : IRange(2, npts + 2 - i - j)) {
                    int l = npts + 3 - i - j - k;
                    double vi = (lob_pts_1d[i - 1] + 1) / 2;
                    double vj = (lob_pts_1d[j - 1] + 1) / 2;
                    double vk = (lob_pts_1d[k - 1] + 1) / 2;
                    double vl = (lob_pts_1d[l - 1] + 1) / 2;


                    localPoints_xi[idx](0) = std::max(0.0, (1 + 3 * vk - vi - vj - vl) / 4);
                    localPoints_xi[idx](1) = std::max(0.0, (1 + 3 * vj - vi - vk - vl) / 4);
                    localPoints_xi[idx](2) = std::max(0.0, (1 + 3 * vi - vj - vk - vl) / 4);

                    ++idx;
                }
            }
        }
        std::sort(localPoints_xi.begin(), localPoints_xi.end(),
                  [](auto const &l, auto const &r) {
                      return std::lexicographical_compare(l.data.rbegin(), l.data.rend(),
                                                          r.data.rbegin(), r.data.rend());
                  });


        int offset = 0;
        cellNodes.reserve(localPoints_xi.size());
        cellTypes.resize(localPoints_xi.size(), CellType::Point1);
        offsets.reserve(localPoints_xi.size() + 1);
        for (auto i : IRange(0, static_cast<int>(localPoints_xi.size()))) {
            cellNodes.push_back(i);
            offsets.push_back(offset++);
        }
        offsets.push_back(offset);

        /*
        int offset = 0;
        idx = 0;
        for (auto layer : IRange(0, npts)) {
            for (auto row : IRange(0, npts - layer)) {
                for (auto i : IRange(0, npts - layer - row)) {

                    //make cell
                    if (layer < npts - 3
                        && row < npts - layer - 3
                        && i < npts - layer - row - 3) {
                        //hex cell
                        int row_stride = npts - layer - row;
                        int layer_stride = (npts - layer) * (npts + 1 - layer) / 2 - row;

                        int idx0 = idx;
                        int idx1 = idx + 1;
                        int idx2 = idx1 + row_stride;
                        int idx3 = idx2 - 1;
                        int idx4 = idx + layer_stride;
                        int idx5 = idx4 + 1;
                        int idx6 = idx5 + row_stride - 1;
                        int idx7 = idx6 - 1;

                        cellNodes.reserve(cellNodes.size() + 8);
                        cellNodes.push_back(idx0);
                        cellNodes.push_back(idx1);
                        cellNodes.push_back(idx2);
                        cellNodes.push_back(idx3);
                        cellNodes.push_back(idx4);
                        cellNodes.push_back(idx5);
                        cellNodes.push_back(idx6);
                        cellNodes.push_back(idx7);

                        cellTypes.push_back(CellType::Hex8);
                        offsets.push_back(offset);
                        offset += 8;
                    } else if (layer < npts - 2
                               && row < npts - layer - 2
                               && i < npts - layer - row - 2) {
                        //make square-base,triangle top tets
                        int row_stride = npts - layer - row;
                        int layer_stride = (npts - layer) * (npts + 1 - layer) / 2 - row;

                        int idx0 = idx;
                        int idx1 = idx + 1;
                        int idx2 = idx1 + row_stride;
                        int idx3 = idx2 - 1;
                        int idx4 = idx + layer_stride;
                        int idx5 = idx4 + 1;
                        int idx6 = idx4 + row_stride - 1;

                        cellNodes.reserve(cellNodes.size() + 4 * 4);

                        cellNodes.push_back(idx1);
                        cellNodes.push_back(idx2);
                        cellNodes.push_back(idx0);
                        cellNodes.push_back(idx5);

                        cellNodes.push_back(idx5);
                        cellNodes.push_back(idx2);
                        cellNodes.push_back(idx0);
                        cellNodes.push_back(idx6);

                        cellNodes.push_back(idx4);
                        cellNodes.push_back(idx6);
                        cellNodes.push_back(idx5);
                        cellNodes.push_back(idx0);

                        cellNodes.push_back(idx3);
                        cellNodes.push_back(idx0);
                        cellNodes.push_back(idx2);
                        cellNodes.push_back(idx6);

                        offsets.push_back(offset);
                        offset += 4;
                        offsets.push_back(offset);
                        offset += 4;
                        offsets.push_back(offset);
                        offset += 4;
                        offsets.push_back(offset);
                        offset += 4;

                        cellTypes.push_back(CellType::Tet4);
                        cellTypes.push_back(CellType::Tet4);
                        cellTypes.push_back(CellType::Tet4);
                        cellTypes.push_back(CellType::Tet4);
                    }
                    else if(layer < npts - 1
                            && row < npts - layer - 1
                            && i < npts - layer - row - 1) {
                        int row_stride = npts - layer - row;
                        int layer_stride = (npts - layer) * (npts + 1 - layer) / 2 - row;

                        int idx0 = idx;
                        int idx1 = idx + 1;
                        int idx2 = idx + row_stride;
                        int idx3 = idx + layer_stride;

                        cellNodes.reserve(cellNodes.size() + 4);
                        cellNodes.push_back(idx0);
                        cellNodes.push_back(idx1);
                        cellNodes.push_back(idx2);
                        cellNodes.push_back(idx3);
                        cellTypes.push_back(CellType::Tet4);
                        offsets.push_back(offset);
                        offset += 4;
                    }


                    //increment idx
                    ++idx;
                }

            }
        }
        offsets.push_back(offset);
        */

    } else {
        //unsupported topoDim
        std::cerr << "Topological Dimension " << elementType.topoDim << std::endl;
    }


    //snap points to help with sorting
    for(auto &x : localPoints_xi) {
        x = snap(x);
    }

    if (elementType.topoDim == 2) {
        quadratureRule.get_triangle_quadrature(std::max(1, 2 * elementType.polyOrder));
    } else if (elementType.topoDim == 3) {
        quadratureRule.get_tetrahedron_quadrature(std::max(1, 2 * elementType.polyOrder));
    }
    localMesh.setGeometryNodes(std::move(localPoints_xi));
    localMesh.setCellNodes(std::move(cellNodes));
    localMesh.setOffsets(std::move(offsets));
    localMesh.setCellTypes(std::move(cellTypes));

    if (elementType.topoDim == 2) {
        triangle_shape_functions(localMesh.getGeometryNodes(),
                                 quadratureRule.nodes,
                                 elementType.polyOrder,
                                 shapeValues,
                                 shapeGradXi);

        build_tri_faces();
    } else if (elementType.topoDim == 3) {
        tetrahedron_shape_functions(localMesh.getGeometryNodes(),
                                    quadratureRule.nodes,
                                    elementType.polyOrder,
                                    shapeValues,
                                    shapeGradXi);
        //build_tet_faces();
    }
}


YAFEL_NAMESPACE_CLOSE