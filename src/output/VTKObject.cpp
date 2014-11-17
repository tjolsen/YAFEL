#include "output/VTKObject.hpp"

YAFEL_NAMESPACE_OPEN

VTKObject::VTKObject(VTKDataType dt, VTKObjectType ot, const std::string &s) : 
  DType(dt), OType(ot), name(s) {}

YAFEL_NAMESPACE_CLOSE
