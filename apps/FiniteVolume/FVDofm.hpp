//
// Created by tyler on 10/15/17.
//

#ifndef YAFEL_FVDOFM_HPP
#define YAFEL_FVDOFM_HPP

#include "yafel_globals.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

/**
 * \class FVDofm
 * DoFManager-like class for finite volume methods
 *
 * Holds the required geometry information for Finite Volume simulations.
 * - Linear Mesh Nodes
 * - List of "full-dimension" (topodim == spatialDim) cells
 *     - Maps FV-indexes to Original element indices
 *     - Also store inverse mapping (with -1 in the invalid element slots)
 * - Volumes of full-dimension elements
 * - Centroids of FV cells
 * - DoFM reference with faces built
 *
 */
class FVDofm {

public:
    FVDofm(DoFManager &dofm_, int8_t NSD);

    DoFManager dofm;
    int8_t nsd;

    std::vector<int> original_cell_index; //map current index -> old_index
    std::vector<int> reverse_index_map;   //map old_index -> current index
    std::vector<coordinate<>> centroids;
    std::vector<double> volumes;
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_FVDOFM_HPP
