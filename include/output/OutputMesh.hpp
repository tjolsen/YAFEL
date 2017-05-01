//
// Created by tyler on 4/11/17.
//

#ifndef YAFEL_OUTPUTMESH_HPP
#define YAFEL_OUTPUTMESH_HPP

#include "yafel_globals.hpp"
#include "utils/DoFManager.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class OutputMesh
{
public:
    OutputMesh(const DoFManager &dm);

    const DoFManager *dofm;


    std::vector<int> local_cells_per_cell;
    std::vector<int> expanded_cells;
    std::vector<int> expanded_cell_offsets;
    std::vector<ElementType> expanded_cell_element_type;
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_OUTPUTMESH_HPP
