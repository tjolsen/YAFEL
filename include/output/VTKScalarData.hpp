#ifndef _YAFEL_VTKSCALARDATA_HPP
#define _YAFEL_VTKSCALARDATA_HPP

#include "yafel_globals.hpp"
#include "output/VTKObject.hpp"
#include "old_handmade_linalg/Vector.hpp"
#include <vector>
#include <string>

YAFEL_NAMESPACE_OPEN

class VTKScalarData : public VTKObject {

private:
  std::vector<double> data;

public:
  VTKScalarData(const std::vector<double> &d, VTKObject::VTKObjectType ot, const std::string &s);
  VTKScalarData(const Vector<double> &d, VTKObject::VTKObjectType ot, const std::string &s);

  void write(FILE *fp);
};

YAFEL_NAMESPACE_CLOSE

#endif
