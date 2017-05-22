//
// Created by tyler on 4/17/17.
//

#include "mesh/CellFace.hpp"
#include "utils/Range.hpp"

YAFEL_NAMESPACE_OPEN

CellFace CellFace::canonicalCellFace(CellType ct, int localFnum)
{
    CellFace ret;
    switch (ct) {
        case CellType::Tri3:
            ret.n_nodes = 2;
            ret.nodes = {{localFnum, (localFnum + 1) % 3}};
            break;
        case CellType::Quad4:
            ret.n_nodes = 2;
            ret.nodes = {{localFnum, (localFnum + 1) % 4}};
            break;

        case CellType::Tet4:
            ret.n_nodes = 3;
            switch (localFnum) {
                case 0:
                    ret.nodes = {{1, 2, 3}};
                    break;
                case 1:
                    ret.nodes = {{0, 3, 2}};
                    break;
                case 2:
                    ret.nodes = {{0, 1, 3}};
                    break;
                case 3:
                    ret.nodes = {{0, 2, 1}};
                    break;
                default:
                    break; //bad news bears if you got here
            }
            break;

        case CellType::Hex8:
            ret.n_nodes = 4;
            switch (localFnum) {
                case 0:
                    ret.nodes = {{0, 4, 7, 3}};
                    break;
                case 1:
                    ret.nodes = {{1, 2, 6, 5}};
                    break;
                case 2:
                    ret.nodes = {{0, 1, 5, 4}};
                    break;
                case 3:
                    ret.nodes = {{2, 3, 7, 6}};
                    break;
                case 4:
                    ret.nodes = {{0, 3, 2, 1}};
                    break;
                case 5:
                    ret.nodes = {{4, 5, 6, 7}};
                    break;
                default:
                    break; //bad news...
            }
        default:
            break;
    }

    return ret;
}

void CellFace::orient()
{
    if (n_nodes == 2) {
        if (nodes[1] < nodes[0]) {
            std::swap(nodes[0], nodes[1]);
            std::swap(left, right);
            std::swap(left_flocal, right_flocal);
            std::swap(left_rot, right_rot);
        }
    } else {
        //n_nodes==3 or n_nodes==4
        int min_i = 0;
        for (auto i : IRange(1, n_nodes)) {
            min_i = (nodes[min_i] < nodes[i]) ? min_i : i;
        }

        left_rot = min_i;

        std::array<int, 4> tmp_nodes(nodes);
        for (auto i : IRange(0, n_nodes)) {
            nodes[i] = tmp_nodes[(min_i + i) % n_nodes];
        }

        if (nodes[n_nodes - 1] < nodes[1]) {
            std::swap(nodes[1], nodes[n_nodes - 1]);
            std::swap(left, right);
            std::swap(left_flocal, right_flocal);
            std::swap(left_rot, right_rot);
        }
    }
}

YAFEL_NAMESPACE_CLOSE