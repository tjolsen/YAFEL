#ifndef _YAFEL_VTKVECTORDATA_HPP
#define _YAFEL_VTKVECTORDATA_HPP

#include "yafel_globals.hpp"
#include "output/VTKObject.hpp"
#include "lin_alg/Vector.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class VTKVectorData : public VTKObject {

private:
  std::vector<Vector> data;

public:
  VTKVectorData(const std::vector<Vector> &d, VTKObject::VTKObjectType ot, const std::string &s);

  void write(FILE *fp);
};

YAFEL_NAMESPACE_CLOSE

#endif
