#ifndef _YAFEL_VTKOUTPUT_HPP
#define _YAFEL_VTKOUTPUT_HPP

#include "yafel_globals.hpp"
#include "output/VTKObject.hpp"
#include "output/VTKMesh.hpp"
#include <vector>
#include <string>

YAFEL_NAMESPACE_OPEN

class VTKOutput {

private:
  VTKObject *vtkmesh;
  std::vector<VTKObject*> pointData;
  std::vector<VTKObject*> cellData;

public:
  VTKOutput();
  
  void addVTKObject(VTKObject* VO);
  void clearData();
  void clearMesh();

  void write(const std::string & fname);
};


YAFEL_NAMESPACE_CLOSE
#endif
