#ifndef _MESHGENERATORHPP
#define _MESHGENERATORHPP

#include "yafel_globals.hpp"
#include "mesh/Mesh.hpp"
#include "lin_alg/Vector.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class MeshGenerator {
  
public:
  using size_type = typename Mesh::size_type;
  using coordinate_type = typename Mesh::coordinate_type;

  MeshGenerator();
  MeshGenerator(size_type NSD);
  MeshGenerator(size_type NSD, double l);
  MeshGenerator(size_type NSD, double l, size_type nx);
  MeshGenerator(size_type NSD, std::vector<double> l);
  MeshGenerator(size_type NSD, std::vector<double> l, std::vector<size_type> nx);
  
  Mesh getMesh();

private:
  size_type nsd;
  std::vector<double> L;
  std::vector<size_type> Nx;

};

YAFEL_NAMESPACE_CLOSE

#endif
