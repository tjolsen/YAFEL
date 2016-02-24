#ifndef __YAFEL_ELEMENTVTKTYPE_HPP
#define __YAFEL_ELEMENTVTKTYPE_HPP

/*
 * Utility (enum) class to provide descriptive names for VTK output
 * element types. The file also includes a mapping from 
 * ElementType to ElementVtkType, which is useful in output operations
 * and eliminates the need for a VTKOutput object to need an ElementFactory,
 * rather than just a Mesh.
 *
 * Names and values for this class are taken from:
 * http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */

#include "yafel_globals.hpp"
#include "utils/ElementType.hpp"


YAFEL_NAMESPACE_OPEN

enum class ElementVtkType : int {
  //linear elements
  VTK_VERTEX=1,
    VTK_POLY_VERTEX=2,
    VTK_LINE=3,
    VTK_POLY_LINE=4,
    VTK_TRIANGLE=5,
    VTK_TRIANGLE_STRIP=6,
    VTK_POLYGON=7,
    VTK_PIXEL=8,
    VTK_QUAD=9,
    VTK_TETRA=10,
    VTK_VOXEL=11,
    VTK_HEXAHEDRON=12,
    VTK_WEDGE=13,
    VTK_PYRAMID=14,
    
  //quadratic elements
    VTK_QUADRATIC_EDGE=21,
    VTK_QUADRATIC_TRIANGLE=22,
    VTK_QUADRATIC_QUAD=23,
    VTK_QUADRATIC_TETRA=24,
    VTK_QUADRATIC_HEXAHEDRON=25,

  //catch-all error type
    VTK_ERROR=-1
};

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

#endif
