#ifndef _MESHGENERATORHPP
#define _MESHGENERATORHPP

#include "yafel_globals.hpp"
#include "mesh/Mesh.hpp"
#include "lin_alg/Vector.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class MeshGenerator {
  
private:
  unsigned nsd;
  std::vector<double> L;
  std::vector<unsigned> Nx;

public:
  MeshGenerator();
  MeshGenerator(unsigned NSD);
  MeshGenerator(unsigned NSD, double l);
  MeshGenerator(unsigned NSD, double l, unsigned nx);
  MeshGenerator(unsigned NSD, std::vector<double> l);
  MeshGenerator(unsigned NSD, std::vector<double> l, std::vector<unsigned> nx);
  
  Mesh getMesh();
};

YAFEL_NAMESPACE_CLOSE

#endif
