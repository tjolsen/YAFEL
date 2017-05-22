//
// Created by tyler on 4/26/17.
//

#include "element/Element.hpp"
#include "element/ShapeFunctionUtils.hpp"


YAFEL_NAMESPACE_OPEN

static std::vector<coordinate<>> make_1d_points(const std::vector<double> &lob_pts_1d);
static std::vector<coordinate<>> make_2d_points(const std::vector<double> &lob_points_1d);
static std::vector<coordinate<>> make_3d_points(const std::vector<double> &lob_points_1d);



void Element::make_simplex()
{
    auto qr = QuadratureRule(elementType.polyOrder + 1,
                             QuadratureRule::QuadratureType::GAUSS_LOBATTO);
    double snap_value = 1.0 / (2 << 20);
    auto snap = [snap_value](auto x) -> coordinate<> { return round(x / snap_value) * snap_value; };

    std::vector<double> lob_pts_1d;
    for (auto x : qr.nodes)
        lob_pts_1d.push_back(x(0));

    std::vector<coordinate<>> localPoints_xi;
    std::vector<int> cellNodes;
    std::vector<int> offsets;
    std::vector<CellType> cellTypes;


    if (elementType.topoDim == 1) {
        //1D simplex is a line, which is also a 1D tensor-product element
        make_tensorProduct();
        return;
    } else if (elementType.topoDim == 2) {

        localPoints_xi = make_2d_points(lob_pts_1d);
        boundaryNodes = make_1d_points(lob_pts_1d);

        //make local triangles
        int nt_local = elementType.polyOrder * elementType.polyOrder;
        cellNodes.resize(3 * nt_local);
        offsets.reserve(nt_local + 1);
        cellTypes.resize(nt_local, CellType::Tri3);

        int idx = 0;
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

        localPoints_xi = make_3d_points(lob_pts_1d);
        boundaryNodes = make_2d_points(lob_pts_1d);

        int offset = 0;
        cellNodes.reserve(localPoints_xi.size());
        cellTypes.resize(localPoints_xi.size(), CellType::Point1);
        offsets.reserve(localPoints_xi.size() + 1);
        for (auto i : IRange(0, static_cast<int>(localPoints_xi.size()))) {
            cellNodes.push_back(i);
            offsets.push_back(offset++);
        }
        offsets.push_back(offset);

    } else {
        //unsupported topoDim
        std::cerr << "Topological Dimension " << elementType.topoDim << std::endl;
    }


    //snap points to help with sorting
    for (auto &x : localPoints_xi) {
        x = snap(x);
    }

    if (elementType.topoDim == 2) {
        quadratureRule.get_triangle_quadrature(std::max(1, 2 * elementType.polyOrder));
        boundaryQuadratureRule = QuadratureRule::make_tensor_product(
                QuadratureRule::QuadratureType::GAUSS_LEGENDRE,
                elementType.topoDim-1,
                2 * elementType.polyOrder);

    } else if (elementType.topoDim == 3) {
        quadratureRule.get_tetrahedron_quadrature(std::max(1, 2 * elementType.polyOrder));
        boundaryQuadratureRule.get_triangle_quadrature(std::max(1, 2 * elementType.polyOrder));
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

        tensor_product_shape_functions(boundaryNodes,
                                       boundaryQuadratureRule.nodes,
                                       elementType.topoDim - 1,
                                       boundaryShapeValues,
                                       boundaryShapeGradXi);

    } else if (elementType.topoDim == 3) {
        tetrahedron_shape_functions(localMesh.getGeometryNodes(),
                                    quadratureRule.nodes,
                                    elementType.polyOrder,
                                    shapeValues,
                                    shapeGradXi);

        triangle_shape_functions(localMesh.getGeometryNodes(),
                                 boundaryQuadratureRule.nodes,
                                 elementType.polyOrder,
                                 boundaryShapeValues,
                                 boundaryShapeGradXi);

    }
}


std::vector<coordinate<>> make_1d_points(const std::vector<double> &lob_pts_1d)
{

    std::vector<coordinate<>> localPoints_xi;
    localPoints_xi.reserve(lob_pts_1d.size());
    for(auto x : lob_pts_1d) {
        localPoints_xi.push_back({x,0,0});
    }
    return localPoints_xi;
}

std::vector<coordinate<>> make_2d_points(const std::vector<double> &lob_pts_1d)
{
    int npts = static_cast<int>(lob_pts_1d.size());
    std::vector<coordinate<>> localPoints_xi;
    // based on scheme be Blyth & Pozrikdis (JAM 2005)
    localPoints_xi.resize((npts * (npts + 1)) / 2);

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
    return localPoints_xi;
}


std::vector<coordinate<>> make_3d_points(const std::vector<double> &lob_pts_1d)
{
    // based on scheme in Blyth & Pozrikdis (JAM 2005)
    int npts = static_cast<int>(lob_pts_1d.size());
    std::vector<coordinate<>> localPoints_xi;
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

    return localPoints_xi;
}


YAFEL_NAMESPACE_CLOSE