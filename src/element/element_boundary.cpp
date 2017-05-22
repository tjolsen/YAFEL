//
// Created by tyler on 5/21/17.
//

#include "element/Element.hpp"
#include <functional>

YAFEL_NAMESPACE_OPEN

std::vector<int> match_face_nodes(const Mesh &elementMesh,
                                  const std::vector<coordinate<>> &boundary_face_coords,
                                  std::function<bool(coordinate<>)> isBoundary,
                                  std::function<coordinate<>(coordinate<>)> faceCoordMapping)
{

    std::vector<int> boundary_nodes;
    std::vector<int> face_nodes(boundary_face_coords.size(), -1);


    auto &elem_node_coords = elementMesh.getGeometryNodes();

    for (int i = 0; i < elementMesh.nNodes(); ++i) {
        if (isBoundary(elem_node_coords[i])) {
            boundary_nodes.push_back(i);
        }
    }


    double MATCH_TOL = 1.0e-6;

    for (auto n : boundary_nodes) {
        auto face_pos = faceCoordMapping(elem_node_coords[n]);

        for (int i = 0; i < face_nodes.size(); ++i) {
            if (norm(face_pos - boundary_face_coords[i]) < MATCH_TOL) {
                face_nodes[i] = n;
            }
        }
    }


    return face_nodes;
}


void Element::build_element_faces()
{

    std::vector<std::function<bool(coordinate<>)>> boundary_functions;
    std::vector<std::function<coordinate<>(coordinate<>)>> boundary_mappings;


    if (elementType.elementTopology == ElementTopology::Simplex && elementType.topoDim == 2) {
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(1)) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(0) + x(1) - 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(0)) < 1.0e-6; });

        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{2*x(0) - 1, 0, 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{2*x(1) - 1, 0, 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{2*(1 - x(1)) - 1, 0}; });

    } else if (elementType.elementTopology == ElementTopology::Simplex && elementType.topoDim == 3) {

        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(0) + x(1) + x(2) - 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(0)) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(1)) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(2)) < 1.0e-6; });

        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(1), x(2), 0.}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(2), x(1), 0.}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(0), x(2), 0.}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(1), x(0), 0.}; });

    } else if (elementType.elementTopology == ElementTopology::TensorProduct && elementType.topoDim == 2) {

        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(1) + 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(0) - 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(1) - 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(0) + 1) < 1.0e-6; });


        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(0), 0, 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(1), 0, 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{-x(0), 0, 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{-x(1), 0, 0}; });

    } else if (elementType.elementTopology == ElementTopology::TensorProduct && elementType.topoDim == 3) {

        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(0) + 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(0) - 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(1) + 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(1) - 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(2) + 1) < 1.0e-6; });
        boundary_functions.emplace_back([](coordinate<> x) { return std::abs(x(2) - 1) < 1.0e-6; });


        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(2), x(1), 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(1), x(2), 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(0), x(2), 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(2), x(0), 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(1), x(0), 0}; });
        boundary_mappings.emplace_back([](coordinate<> x) { return coordinate<>{x(0), x(1), 0}; });
    }

    for (int face = 0; face < boundary_mappings.size(); ++face) {

        auto tmp = match_face_nodes(
                localMesh,
                boundaryNodes,
                boundary_functions[face],
                boundary_mappings[face]
        );
        face_nodes.push_back(tmp);

    }
}


YAFEL_NAMESPACE_CLOSE