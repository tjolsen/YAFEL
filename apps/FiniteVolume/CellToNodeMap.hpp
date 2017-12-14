//
// Created by tyler on 10/15/17.
//

#ifndef YAFEL_INTERPOLATECELLSTONODES_HPP
#define YAFEL_INTERPOLATECELLSTONODES_HPP

#include "yafel_globals.hpp"
#include "FVDofm.hpp"

#include <Eigen/Sparse>
#include <Eigen/Core>

YAFEL_NAMESPACE_OPEN

template<int NSD>
class CellToNodeMap {

public:
    CellToNodeMap(FVDofm const& fvdofm) : fvDofm(fvdofm) { construct(); }

    void interpolateToNodes(Eigen::VectorXd const& src, Eigen::VectorXd &dst) {
        if(!is_constructed) {
            construct();
        }
        dst.noalias() = interpolation_matrix*src;
    }

private:
    const FVDofm & fvDofm;
    Eigen::SparseMatrix<double,Eigen::RowMajor> interpolation_matrix;
    void construct();
    bool is_constructed{false};
};



template<int NSD>
void CellToNodeMap<NSD>::construct()
{
    std::vector<double> node_adjacent_inverse_volumes(fvDofm.dofm.nNodes(),0.0);
    std::vector<Eigen::Triplet<double>> adjacentTriplets;

    std::vector<int> element_container;
    for(int cell=0; cell < fvDofm.centroids.size(); ++cell) {
        auto cell_idx = fvDofm.original_cell_index[cell];
        auto v = fvDofm.volumes[cell_idx];
        auto Vinv_cell = 1.0/(v == 0.0 ? -1.0 : v);
        fvDofm.dofm.getGlobalNodes(cell_idx, element_container);
        for(auto n : element_container) {
            auto& triplet = adjacentTriplets.emplace_back(n, cell_idx, Vinv_cell);
            node_adjacent_inverse_volumes[n] += Vinv_cell;
        }
    }

    for(auto& t : adjacentTriplets) {
        t = Eigen::Triplet<double>{t.row(), t.col(), t.value()/node_adjacent_inverse_volumes[t.row()]};
    }

    interpolation_matrix.resize(node_adjacent_inverse_volumes.size(), fvDofm.volumes.size());
    interpolation_matrix.setFromTriplets(adjacentTriplets.begin(), adjacentTriplets.end());
}

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_INTERPOLATECELLSTONODES_HPP
