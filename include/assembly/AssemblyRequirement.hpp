//
// Created by tyler on 5/2/17.
//

#ifndef YAFEL_ASSEMBLYREQUIREMENT_HPP
#define YAFEL_ASSEMBLYREQUIREMENT_HPP

#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

enum class AssemblyRequirement : int
{
    Residual,
    Tangent,
    DtMass,
    DtDtMass
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_ASSEMBLYREQUIREMENT_HPP
