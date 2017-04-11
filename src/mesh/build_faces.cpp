//
// Created by tyler on 4/6/17.
//

#include "mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include <algorithm>

YAFEL_NAMESPACE_OPEN

void Mesh::buildInternalFaces()
{

    //build naive list of faces
    int num_naive_faces{0};
    for (auto i : IRange(0, nCells())) {
        int n{0};
        switch (getCellType(i)) {
            case CellType::Tri3:
                n = 3;
                break;
            case CellType::Tet4:
            case CellType::Quad4:
                n = 4;
                break;
            case CellType::Hex8:
                n = 6;
                break;
            default:
                n = 0;
        }

        num_naive_faces += n;
    }
    internal_faces_.clear();
    internal_faces_.reserve(num_naive_faces);
    std::vector<int> cell_nodes(8);

    for (auto i : IRange(0, nCells())) {

        getCellNodes(i, cell_nodes);

        switch (getCellType(i)) {
            case CellType::Tri3:
                for (auto f : IRange(0, 3)) {
                    internal_faces_.emplace_back();
                    auto &F = internal_faces_.back();
                    F.n_nodes = 2;
                    F.left = i;
                    F.right = -1;
                    F.nodes = {cell_nodes[f], cell_nodes[(f + 1) % 3]};
                    F.left_flocal = f;
                    F.orient();
                }
                break;
            case CellType::Quad4:
                for (auto f : IRange(0, 4)) {
                    internal_faces_.emplace_back();
                    auto &F = internal_faces_.back();
                    F.n_nodes = 2;
                    F.left = i;
                    F.right = -1;
                    F.nodes = {cell_nodes[f], cell_nodes[(f + 1) % 4]};
                    F.left_flocal = f;
                    F.orient();
                }
                break;
            default:
                break;
        }

    }

    std::sort(internal_faces_.begin(), internal_faces_.end(),
              [](const CellFace &L, const CellFace &R) { return L.nodes < R.nodes; }
    );

    for (auto i : IRange(1, static_cast<int>(internal_faces_.size()))) {
        auto &F = internal_faces_[i];
        auto &Fprev = internal_faces_[i - 1];
        if (F.nodes == Fprev.nodes) {

            Fprev.left = std::max(Fprev.left, F.left);
            Fprev.right = std::max(Fprev.right, F.right);
            Fprev.left_flocal = std::max(Fprev.left_flocal, F.left_flocal);
            Fprev.right_flocal = std::max(Fprev.right_flocal, F.right_flocal);

        }

    }

    auto new_end = std::unique(internal_faces_.begin(),
                               internal_faces_.end(),
                               [](const CellFace &L, const CellFace &R) { return L.nodes == R.nodes; });

    auto new_size = std::distance(internal_faces_.begin(), new_end);
    internal_faces_.resize(new_size);
    internal_faces_.shrink_to_fit();
}

YAFEL_NAMESPACE_CLOSE