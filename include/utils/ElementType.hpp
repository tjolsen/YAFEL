#ifndef _YAFEL_ELEMENTTYPE_HPP
#define _YAFEL_ELEMENTTYPE_HPP

/*
 * Utility (enum) class to provide a unified way to describe an element
 * type. This file also provides a set of functions to map from a
 * preprocessor element type description to the one used here.
 */


#include "yafel_globals.hpp"
#include <cstddef>

YAFEL_NAMESPACE_OPEN

enum class ElementType {
  // 1-D Elements
  LINEAR_LINE,
    QUADRATIC_LINE,
    CUBIC_LINE,
    

  // 2-D Quadrilateral Elements
    LINEAR_QUAD,
    QUADRATIC_QUAD,
    CUBIC_QUAD,
    
  // 2-D Triangle Elements
    LINEAR_TRI,
    QUADRATIC_TRI,
    CUBIC_TRI,
    
  // 3-D Hexahedral Elements
    LINEAR_HEX,
    QUADRATIC_HEX,
    CUBIC_HEX,
    
  // 3-D Tetrahedral Elements
    LINEAR_TET,
    QUADRATIC_TET,
    CUBIC_TET,

  // DG Types
    DG_QUAD,
    
  // Error Type
    NULL_ELEMENT
};

//functions for mapping preprocessor-specific element types to yafel representation
namespace ElementType_Mappings{
  ElementType gmsh_to_ElementType(std::size_t eltype);
}  

YAFEL_NAMESPACE_CLOSE


#endif
