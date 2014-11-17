#ifndef _YAFEL_VTKOUTPUT_HPP
#define _YAFEL_VTKOUTPUT_HPP

#include <vector>
#include "yafel_globals.hpp"
#include "output/VTKObject.hpp"
#include "output/VTKMesh.hpp"

YAFEL_NAMESPACE_OPEN

class VTKOutput {

private:
  VTKObject *vtkmesh;
  std::vector<VTKObject*> pointData;
  std::vector<VTKObject*> cellData;

public:
  VTKOutput();
  void addVTKObject(VTKObject * VO);
  void write(const char* fname);
};

YAFEL_NAMESPACE_CLOSE
#endif
