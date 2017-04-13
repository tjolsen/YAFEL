//
// Created by tyler on 4/5/17.
//

#ifndef YAFEL_CELLFACE_HPP
#define YAFEL_CELLFACE_HPP

#include "yafel_globals.hpp"
#include "utils/Range.hpp"
#include <array>
#include <iostream>

YAFEL_NAMESPACE_OPEN

/**
 * \class CellFace
 *
 * \brief Structure to represent a face between cells
 *
 * Intended to be used to represent boundaries between all types of cells.
 * The "nodes" container can hold up to 4 node indices.
 * Un-set values for left, right are defined to be -1.
 * This will be useful for ultimately determining whether a face
 * comprises a global boundary.
 *
 * For boundaries between 2D elemnets (quad/tri), the "nodes"
 * container will hold 2 integers corresponding to the corners between
 * the cells.
 * The orientation of the face is defined so that nodes[0] < nodes[1].
 *
 * For 3D elements (tet/hex), the "nodes" container will
 * hold 3/4 nodes (as appropriate), and the orientation is
 * defined so that nodes[0] = min(nodes), nodes[1] < nodes[2].
 * This uniquely defines the face orientation regardless
 * of how it occurs in the original mesh.
 */
struct CellFace {

    inline CellFace()
            : nodes(), left(-1), right(-1), n_nodes(0),
              left_flocal(-1), right_flocal(-1)
    {
        this->nodes.fill(0);
    }

    //node idx container
    std::array<int,4> nodes;

    //Indices of left/right cells
    int left;
    int right;
    int n_nodes;

    //Local face indices
    int left_flocal;
    int right_flocal;

    inline void orient() {
        if(n_nodes == 2) {
            if(nodes[1] < nodes[0]) {
                std::swap(nodes[0], nodes[1]);
                std::swap(left,right);
                std::swap(left_flocal,right_flocal);
            }
        }
        else {
            //n_nodes==3 or n_nodes==4
            int min_i=0;
            for(auto i : IRange(1,n_nodes)) {
                min_i = (nodes[min_i] < nodes[i]) ? min_i : i;
            }

            std::array<int,4> tmp_nodes(nodes);
            for(auto i : IRange(0,n_nodes)) {
                nodes[i] = tmp_nodes[(min_i+i)%n_nodes];
            }

            if(nodes[n_nodes-1] < nodes[1]) {
                std::swap(nodes[1], nodes[n_nodes-1]);
                std::swap(left,right);
                std::swap(left_flocal,right_flocal);
            }
        }
    };
};

inline std::ostream& operator<<(std::ostream& out, const CellFace& F)
{
    out << "F{ {";
    for(int i=0; i<F.n_nodes-1; ++i) {
        out << F.nodes[i] << ", ";
    }
    out << F.nodes[F.n_nodes-1]
        << "}, L(" << F.left
        << "), R(" << F.right
        << "), Lf(" << F.left_flocal
        << "), Rf(" << F.right_flocal
        << ") }";

    return out;
}

YAFEL_NAMESPACE_CLOSE


#endif //YAFEL_CELLFACE_HPP
