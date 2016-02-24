#ifndef _YAFEL_MESHTOPOLOGY_HPP
#define _YAFEL_MESHTOPOLOGY_HPP

#include "mesh/GenericMesh.hpp"
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
  std::map<std::size_t, TopoFace*> faces;

  std::size_t get_line_id(std::size_t vid1, std::size_t vid2, bool &exists);

public:
  template<typename MT, unsigned NSD>
  MeshTopology(const GenericMesh<MT,NSD> &M);
  ~MeshTopology();
  void print(std::ostream &out);

  // method to return the id of the edgenum'th neighbor of cell facenum.
  // returns facenum on any failure
  std::size_t getCellNeighbor(std::size_t faceNum, std::size_t edgenum);
};

YAFEL_NAMESPACE_CLOSE

#endif
