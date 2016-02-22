#ifndef __YAFEL_ELEMENTTYPE_HPP
#define __YAFEL_ELEMENTTYPE_HPP

#include "yafel_globals.hpp"

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

  // Error Type
    NULL_ELEMENT
};


namespace ElementType_Mappings{
  
  ElementType gmsh_to_ElementType(std::size_t eltype) {
    switch(eltype) {
    case 1:
      return ElementType::LINEAR_LINE;
    case 2:
      return ElementType::LINEAR_TRI;
    case 3:
      return ElementType::LINEAR_QUAD; 
    case 4:
      return ElementType::LINEAR_TET;
    case 5:
      return ElementType::LINEAR_HEX;
    case 8:
      return ElementType::QUADRATIC_LINE;
    case 9:
      return ElementType::QUADRATIC_TRI;
    case 10:
      return ElementType::QUADRATIC_QUAD;
    case 11:
      return ElementType::QUADRATIC_TET;
    case 12:
      return ElementType::QUADRATIC_HEX;
    default:
      return ElementType::NULL_ELEMENT;
    }
  }

}


YAFEL_NAMESPACE_CLOSE


#endif
