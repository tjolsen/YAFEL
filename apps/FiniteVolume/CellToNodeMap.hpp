//
// Created by tyler on 10/15/17.
//

#ifndef YAFEL_INTERPOLATECELLSTONODES_HPP
#define YAFEL_INTERPOLATECELLSTONODES_HPP

#include "yafel_globals.hpp"
#include "FVDofm.hpp"

#include <eigen3/Eigen/Sparse>

YAFEL_NAMESPACE_OPEN

template<int NSD>
class CellToNodeMap {

public:
    CellToNodeMap(FVDofm const& fvdofm) : fvDofm(fvdofm) { construct(); }

private:
    const FVDofm & fvDofm;
    Eigen::SparseMatrix<double,Eigen::RowMajor> interpolation_matrix;
    void construct();
};



template<int NSD>
void CellToNodeMap<NSD>::construct()
{
    std::vector<double> node_adjacent_inverse_volumes(fvDofm.dofm.nNodes(),0.0);
    std::vector<Eigen::Triplet<double>> adjacentTriplets;

    std::vector<int> element_container;
    for(int cell=0; cell < fvDofm.centroids.size(); ++cell) {
        auto Vinv_cell = 1.0/fvDofm.volumes[cell];
        fvDofm.dofm.getGlobalNodes(fvDofm.original_cell_index[cell], element_container);
        for(auto n : element_container) {
            auto& triplet = adjacentTriplets.emplace_back(cell, n, Vinv_cell);
            node_adjacent_inverse_volumes[n] += Vinv_cell;
        }
    }

    for(auto &t : adjacentTriplets) {
        t.value() /= node_adjacent_inverse_volumes[t.col()];
    }

    interpolation_matrix.resize(fvDofm.volumes.size(), node_adjacent_inverse_volumes.size());
    interpolation_matrix.setFromTriplets(adjacentTriplets.begin(), adjacentTriplets.end());
}

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_INTERPOLATECELLSTONODES_HPP
