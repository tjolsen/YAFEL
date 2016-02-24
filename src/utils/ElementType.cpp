#include "utils/ElementType.hpp"
#include <cstddef>

YAFEL_NAMESPACE_OPEN

namespace ElementType_Mappings {

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
  
  
}//end namespace ElementType_Mappings




YAFEL_NAMESPACE_CLOSE
