//
// Created by tyler on 4/6/17.
//

#include "mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include <algorithm>

YAFEL_NAMESPACE_OPEN

static int numCellFaces(CellType ct) {
    switch (ct) {
        case CellType::Tri3:
            return 3;
        case CellType::Tet4:
        case CellType::Quad4:
            return 4;
        case CellType::Hex8:
            return 6;
        default:
            return 0;
    }

}

void Mesh::buildInternalFaces()
{

    //build naive list of faces
    int num_naive_faces{0};
    for (auto i : IRange(0, nCells())) {
        num_naive_faces += numCellFaces(getCellType(i));
    }

    internal_faces_.clear();
    internal_faces_.reserve(num_naive_faces);
    std::vector<int> cell_nodes(8);

    for (auto i : IRange(0, nCells())) {

        getCellNodes(i, cell_nodes);
        auto ct = getCellType(i);

        for(auto f : IRange(0, numCellFaces(ct))) {
            auto F = CellFace::canonicalCellFace(ct,f);
            for(auto n : IRange(0,F.n_nodes)) {
                F.nodes[n] = cell_nodes[F.nodes[n]];
            }
            F.left = i;
            F.left_flocal = f;
            F.orient();
            internal_faces_.push_back(F);
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

    for (auto i : IRange(0, static_cast<int>(internal_faces_.size()))) {
        if (internal_faces_[i].left < 0 || internal_faces_[i].right < 0) {
            boundary_face_idxs_.push_back(i);
        }
    }
}

YAFEL_NAMESPACE_CLOSE