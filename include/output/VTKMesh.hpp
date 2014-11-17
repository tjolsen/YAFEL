#ifndef _YAFEL_VTKMESH_HPP
#define _YAFEL_VTKMESH_HPP

#include <cstdio>
#include "yafel_globals.hpp"
#include "element/ElementFactory.hpp"
#include "output/VTKObject.hpp"

YAFEL_NAMESPACE_OPEN

class VTKMesh : public VTKObject {
  
private:
  ElementFactory * EFp;

public:
  VTKMesh(ElementFactory &EF);
  void write(FILE *fp);
};

YAFEL_NAMESPACE_CLOSE

#endif
