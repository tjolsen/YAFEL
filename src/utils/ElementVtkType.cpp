#include "yafel_globals.hpp"
#include "utils/ElementVtkType.hpp"

YAFEL_NAMESPACE_OPEN

ElementVtkType ElementType_to_ElementVtkType(ElementType e) {
  
  switch(e) {
  case ElementType::LINEAR_LINE:
    return ElementVtkType::VTK_LINE;

  case ElementType::QUADRATIC_LINE:
    return ElementVtkType::VTK_QUADRATIC_EDGE;

  case ElementType::LINEAR_QUAD:
    return ElementVtkType::VTK_QUAD;

  case ElementType::QUADRATIC_QUAD:
    return ElementVtkType::VTK_QUADRATIC_QUAD;

  case ElementType::LINEAR_TRI:
    return ElementVtkType::VTK_TRIANGLE;

  case ElementType::QUADRATIC_TRI:
    return ElementVtkType::VTK_QUADRATIC_TRIANGLE;

  case ElementType::LINEAR_HEX:
    return ElementVtkType::VTK_HEXAHEDRON;

  case ElementType::QUADRATIC_HEX:
    return ElementVtkType::VTK_QUADRATIC_HEXAHEDRON;

  case ElementType::LINEAR_TET:
    return ElementVtkType::VTK_TETRA;

  case ElementType::QUADRATIC_TET:
    return ElementVtkType::VTK_QUADRATIC_TETRA;

  default:
    return ElementVtkType::VTK_ERROR;
  }
}

YAFEL_NAMESPACE_CLOSE
