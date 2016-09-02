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



std::string to_string(ElementType ET) {
    switch(ET) {
    // 1-D Elements
    case ElementType::LINEAR_LINE:
	return "LINEAR_LINE";
    case ElementType::QUADRATIC_LINE:
	return "QUADRATIC_LINE";
    case ElementType::CUBIC_LINE:
	return "CUBIC_LINE";

    // 2-D Quadrilateral Elements
    case ElementType::LINEAR_QUAD:
	return "LINEAR_QUAD";
    case ElementType::QUADRATIC_QUAD:
	return "QUADRATIC_QUAD";
    case ElementType::CUBIC_QUAD:
	return "CUBIC_QUAD";
    
    // 2-D Triangle Elements
    case ElementType::LINEAR_TRI:
	return "LINEAR_TRI";
    case ElementType::QUADRATIC_TRI:
	return "QUADRATIC_TRI";
    case ElementType::CUBIC_TRI:
	return "CUBIC_TRI";
    
    // 3-D Hexahedral Elements
    case ElementType::LINEAR_HEX:
	return "LINEAR_HEX";
    case ElementType::QUADRATIC_HEX:
	return "QUADRATIC_HEX";
    case ElementType::CUBIC_HEX:
	return "CUBIC_HEX";
    
    // 3-D Tetrahedral Elements
    case ElementType::LINEAR_TET:
	return "LINEAR_TET";
    case ElementType::QUADRATIC_TET:
	return "QUADRATIC_TET";
    case ElementType::CUBIC_TET:
	return "CUBIC_TET";

    // DG Types
    case ElementType::DG_QUAD:
	return "DG_QUAD";
    
    // Error Type
    case ElementType::NULL_ELEMENT:
    default:
	return "NULL_ELEMENT";

    }
}





YAFEL_NAMESPACE_CLOSE
