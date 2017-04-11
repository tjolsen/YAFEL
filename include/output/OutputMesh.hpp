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
    inline OutputMesh(const DoFManager &dm) : dofm(&dm) {}

    const DoFManager *dofm;

    std::vector<int> local_cells_per_cell;
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_OUTPUTMESH_HPP
