#ifndef _YAFEL_VTKTENSORDATA_HPP
#define _YAFEL_VTKTENSORDATA_HPP

#include "yafel_globals.hpp"
#include "output/VTKObject.hpp"
#include "lin_alg/FullMatrix.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

// WARNING, VTK CURRENTLY ONLY SUPPORTS 3x3 SYMMETRIC TENSORS

class VTKTensorData : public VTKObject {

private:
  std::vector<FullMatrix> data;

public:
  VTKTensorData(const std::vector<FullMatrix> &d, VTKObject::VTKObjectType ot, const std::string &s);

  void write(FILE *fp);
};

YAFEL_NAMESPACE_CLOSE

#endif
