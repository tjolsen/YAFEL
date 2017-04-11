//
// Created by tyler on 4/11/17.
//

#ifndef YAFEL_OUTPUTMESH_HPP
#define YAFEL_OUTPUTMESH_HPP

#include "yafel_globals.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

class OutputMesh
{
public:
    inline OutputMesh(const DoFManager &dm) : dofm(&dm) {}

    const DoFManager *dofm;
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_OUTPUTMESH_HPP
