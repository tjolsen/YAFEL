//
// Created by tyler on 4/18/17.
//

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "element/element_boundary_nodes.hpp"

YAFEL_NAMESPACE_OPEN

void Element::build_tet_faces()
{
    std::vector<std::function<bool(coordinate<>)>> boundary_funcs{
            [](coordinate<> xi) { return std::abs(sum(xi) - 1.0) < 1.0e-6; },
            [](coordinate<> xi) { return std::abs(xi(0)) < 1.0e-6; },
            [](coordinate<> xi) { return std::abs(xi(1)) < 1.0e-6; },
            [](coordinate<> xi) { return std::abs(xi(2)) < 1.0e-6; },
    };

    face_perm.resize(4);
    for (auto &F : face_perm) {
        F.resize(3);
        for (auto &f : F) {
            f.resize(2);
        }
    }


    std::vector<std::vector<int>> face_node_indices(4);
    for(auto n_idx : IRange(0,localMesh.nNodes())) {
        auto xi = localMesh.getGeometryNodes()[n_idx];
        for(auto f_idx : IRange(0,4)) {
            if(boundary_funcs[f_idx](xi)) {
                face_node_indices[f_idx].push_back(n_idx);
            }
        }
    }

    std::vector<std::function<coordinate<>(coordinate<>)>> space_to_face_coords
            {
                    [](coordinate<> xi)->coordinate<>{return {xi(1),xi(2),0};},
                    [](coordinate<> xi)->coordinate<>{return {xi(1), xi(2),0};},
                    [](coordinate<> xi)->coordinate<>{return {xi(0), xi(2),0};},
                    [](coordinate<> xi)->coordinate<>{return {xi(0), xi(1), 0};}
            };


    std::vector<std::function<coordinate<>(coordinate<>)>> face_to_neighbor_coords
            {
                    [](coordinate<> xi)->coordinate<>{return {xi(1),xi(0),0};},
                    [](coordinate<> xi)->coordinate<>{return {1 - sum(xi),xi(1),0};},
                    [](coordinate<> xi)->coordinate<>{return {xi(0),1-sum(xi),0};}
            };



    for(auto f_idx : IRange(0,4)) {

        std::vector<coordinate<>> f_coords;
        for(auto n : face_node_indices[f_idx]) {
            auto x = localMesh.getGeometryNodes()[n];
            auto xf = space_to_face_coords[f_idx](x);
            f_coords.push_back(xf);
        }

        for(auto p : IRange(0,3)) {
            std::vector<coordinate<>> neighbor_f_coords;

            for(auto xf : f_coords) {
                auto nxf = face_to_neighbor_coords[p](xf);
                neighbor_f_coords.push_back(nxf);
            }

            //do something with the neighbor face params (sorting?)
        }

    }
}


YAFEL_NAMESPACE_CLOSE