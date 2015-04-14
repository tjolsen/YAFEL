#ifndef _YAFEL_MESHTOPOLOGY_HPP
#define _YAFEL_MESHTOPOLOGY_HPP

#include "mesh/Mesh.hpp"
#include "mesh/TopoPoint.hpp"
#include "mesh/TopoLine.hpp"
#include "mesh/TopoFace.hpp"
//#include "mesh/TopoCell.hpp"

#include <vector>
#include <map>
#include <iostream>

YAFEL_NAMESPACE_OPEN

class MeshTopology {

private:
  std::vector<TopoPoint> points;
  std::vector<TopoLine*> lines;
  std::map<unsigned, TopoFace*> faces;

  unsigned get_line_id(unsigned vid1, unsigned vid2, bool &exists);

public:
  MeshTopology(const Mesh &M);
  ~MeshTopology();
  void print(std::ostream &out);
  
};

YAFEL_NAMESPACE_CLOSE

#endif
