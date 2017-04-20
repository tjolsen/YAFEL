//
// Created by tyler on 4/18/17.
//

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "element/element_boundary_nodes.hpp"

YAFEL_NAMESPACE_OPEN

void Element::build_tet_faces()
{

    Tensor<3, 1> n0{1, 1, 1};
    Tensor<3, 1> n1{-1, 0, 0};
    Tensor<3, 1> n2{0, -1, 0};
    Tensor<3, 1> n3{0, 0, -1};
    std::vector<Tensor<3, 1>> normals{n0, n1, n2, n3};

    Tensor<3, 1> e0_0F{-1, 0, 1}, e0_1F{0, 1, -1}, e0_2F{1, -1, 0};
    Tensor<3, 1> e0_0R(-e0_2F), e0_1R(-e0_0F), e0_2R(-e0_1F);

    Tensor<3, 1> e1_0F{0, 0, 1}, e1_1F{0, 1, -1}, e1_2F{0, -1, 0};
    Tensor<3, 1> e1_0R(-e0_2F), e1_1R(-e0_0F), e1_2R(-e0_1F);

    Tensor<3, 1> e2_0F{1, 0, -1}, e2_1F{-1, 0, 0}, e2_2F{0, 0, 1};
    Tensor<3, 1> e2_0R(-e0_2F), e2_1R(-e0_0F), e2_2R(-e0_1F);

    Tensor<3, 1> e3_0F{1, 0, 0}, e3_1F{-1, 1, 0}, e3_2F{0, -1, 0};
    Tensor<3, 1> e3_0R(-e0_2F), e3_1R(-e0_0F), e3_2R(-e0_1F);


    std::vector<Tensor<3, 1>> e0F{e0_0F, e0_1F, e0_2F};
    std::vector<Tensor<3, 1>> e0R{e0_0R, e0_1R, e0_2R};
    std::vector<Tensor<3, 1>> e1F{e1_0F, e1_1F, e1_2F};
    std::vector<Tensor<3, 1>> e1R{e1_0R, e1_1R, e1_2R};
    std::vector<Tensor<3, 1>> e2F{e2_0F, e2_1F, e2_2F};
    std::vector<Tensor<3, 1>> e2R{e2_0R, e2_1R, e2_2R};
    std::vector<Tensor<3, 1>> e3F{e3_0F, e3_1F, e3_2F};
    std::vector<Tensor<3, 1>> e3R{e3_0R, e3_1R, e3_2R};

    std::vector<std::vector<Tensor<3, 1>>> eF{e0F, e1F, e2F, e3F};
    std::vector<std::vector<Tensor<3, 1>>> eR{e0R, e1R, e2R, e3R};

    std::vector<std::function<bool(coordinate<>)>> boundary_funcs{
            [](coordinate<> xi) { return std::abs(dot(xi, Tensor<3, 1>{1, 1, 1}) - 1.0) < 1.0e-6; },
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

    for (auto face : IRange(0, 4)) {
        for (auto rot : IRange(0, 3)) {

            auto face_nodes = make_element_boundary_nodes(localMesh.getGeometryNodes(),
                                                          normals[face], eF[face][rot],
                                                          boundary_funcs[face]);
            face_perm[face][rot][0].swap(face_nodes);


            face_nodes = make_element_boundary_nodes(localMesh.getGeometryNodes(),
                                                     -normals[face], eF[face][rot],
                                                     boundary_funcs[face]);
            face_perm[face][rot][1].swap(face_nodes);
        }
    }
}


YAFEL_NAMESPACE_CLOSE