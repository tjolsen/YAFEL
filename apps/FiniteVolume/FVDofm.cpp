//
// Created by tyler on 10/15/17.
//

#include "FVDofm.hpp"
#include "cell_centroids.hpp"

YAFEL_NAMESPACE_OPEN

FVDofm::FVDofm(DoFManager &dofm_, int8_t NSD)
        : dofm(dofm_), nsd(NSD)
{

    reverse_index_map.resize(dofm.nCells(), -1);
    int n_FV_cells{0};
    for(auto i=0; i<dofm.nCells(); ++i) {
        if(dofm.element_types[i].topoDim == NSD) {
            ++n_FV_cells;
        }
    }

    original_cell_index.resize(n_FV_cells);
    int idx{0};
    for(int i=0; i<dofm.nCells(); ++i) {
        if(dofm.element_types[i].topoDim == NSD) {
            original_cell_index[idx] = i;
            reverse_index_map[i] = idx;

            ++idx;
        }
    }

    switch(NSD) {
        case 1:
            auto [tmp_centroids, tmp_volumes] = cell_cenetroids<1>(); break;
        case 2:
            auto [tmp_centroids, tmp_volumes] = cell_cenetroids<2>(); break;
        case 3:
            auto [tmp_centroids, tmp_volumes] = cell_cenetroids<3>(); break;
    }
    centroids = std::move(tmp_centroids);
    volumes = std::move(tmp_volumes);

}


YAFEL_NAMESPACE_CLOSE
