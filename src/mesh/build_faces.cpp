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
                    F.nodes = {{cell_nodes[f], cell_nodes[(f + 1) % 3]}};
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
                    F.nodes = {{cell_nodes[f], cell_nodes[(f + 1) % 4]}};
                    F.left_flocal = f;
                    F.orient();
                }
                break;

            case CellType::Tet4: {
                CellFace F0, F1, F2, F3;
                F0.left = i;
                F0.right = -1;
                F0.left_flocal = 0;
                F0.nodes = {{cell_nodes[1], cell_nodes[3], cell_nodes[2]}};
                F0.n_nodes = 3;
                F1.left = i;
                F1.right = -1;
                F1.left_flocal = 1;
                F1.nodes = {{cell_nodes[0], cell_nodes[2], cell_nodes[3]}};
                F1.n_nodes = 3;
                F2.left = i;
                F2.right = -1;
                F2.left_flocal = 2;
                F2.nodes = {{cell_nodes[0], cell_nodes[3], cell_nodes[1]}};
                F2.n_nodes = 3;
                F3.left = i;
                F3.right = -1;
                F3.left_flocal = 3;
                F3.nodes = {{cell_nodes[0], cell_nodes[1], cell_nodes[2]}};
                F3.n_nodes = 3;
                F0.orient();
                F1.orient();
                F2.orient();
                F3.orient();
                internal_faces_.push_back(F0);
                internal_faces_.push_back(F1);
                internal_faces_.push_back(F2);
                internal_faces_.push_back(F3);
            }
                break;

            case CellType::Hex8: {
                {
                    CellFace F;
                    F.n_nodes = 4;
                    F.left = i;
                    F.right = -1;
                    F.n_nodes = 4;
                    F.left_flocal = 0;
                    F.right_flocal = -1;
                    F.nodes = {{cell_nodes[0],
                                       cell_nodes[3],
                                       cell_nodes[7],
                                       cell_nodes[4]}};
                    F.orient();
                    internal_faces_.push_back(F);
                }
                {
                    CellFace F;
                    F.n_nodes = 4;
                    F.left = i;
                    F.right = -1;
                    F.n_nodes = 4;
                    F.left_flocal = 1;
                    F.right_flocal = -1;
                    F.nodes = {{cell_nodes[1],
                                       cell_nodes[5],
                                       cell_nodes[6],
                                       cell_nodes[2]}};
                    F.orient();
                    internal_faces_.push_back(F);
                }

                {
                    CellFace F;
                    F.n_nodes = 4;
                    F.left = i;
                    F.right = -1;
                    F.n_nodes = 4;
                    F.left_flocal = 2;
                    F.right_flocal = -1;
                    F.nodes = {{cell_nodes[0],
                                       cell_nodes[4],
                                       cell_nodes[5],
                                       cell_nodes[1]}};
                    F.orient();
                    internal_faces_.push_back(F);
                }

                {
                    CellFace F;
                    F.n_nodes = 4;
                    F.left = i;
                    F.right = -1;
                    F.left_flocal = 3;
                    F.right_flocal = -1;
                    F.nodes = {{cell_nodes[3],
                                       cell_nodes[2],
                                       cell_nodes[6],
                                       cell_nodes[7]}};
                    F.orient();
                    internal_faces_.push_back(F);
                }

                {
                    CellFace F;
                    F.n_nodes = 4;
                    F.left = i;
                    F.right = -1;
                    F.n_nodes = 4;
                    F.left_flocal = 4;
                    F.right_flocal = -1;
                    F.nodes = {{cell_nodes[0],
                                       cell_nodes[1],
                                       cell_nodes[2],
                                       cell_nodes[3]}};
                    F.orient();
                    internal_faces_.push_back(F);
                }

                {
                    CellFace F;
                    F.n_nodes = 4;
                    F.left = i;
                    F.right = -1;
                    F.n_nodes = 4;
                    F.left_flocal = 5;
                    F.right_flocal = -1;
                    F.nodes = {{cell_nodes[5],
                                       cell_nodes[4],
                                       cell_nodes[7],
                                       cell_nodes[6]}};
                    F.orient();
                    internal_faces_.push_back(F);
                }
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

    for (auto i : IRange(0, static_cast<int>(internal_faces_.size()))) {
        if (internal_faces_[i].left < 0 || internal_faces_[i].right < 0) {
            boundary_face_idxs_.push_back(i);
        }
    }
}

YAFEL_NAMESPACE_CLOSE