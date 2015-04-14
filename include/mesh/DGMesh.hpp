#ifndef _YAFEL_DGMESH_HPP
#define _YAFEL_DGMESH_HPP

#include "yafel_globals.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/MeshTopology.hpp"

YAFEL_NAMESPACE_OPEN

/*
  DGMesh is currently ONLY GOING TO WORK FOR 2-MANIFOLD MESHES!!!
  I need to think about how to handle lines in a meaningful way in a
  volume mesh.
 */


class DGMesh : public Mesh {

private:
  MeshTopology meshTopology;
  
public:
  DGMesh();
  DGMesh(const std::vector<Vector> &nodes,
	 const std::vector<std::vector<unsigned> &elems,
	 const std::vector<unsigned> &eltype);
  DGMesh(const std::vector<Vector> &nodes,
	 const std::vector<std::vector<unsigned> > &elems,
	 const std::vector<unsigned> &eltype,
	 const std::vector<std::vector<unsigned> > &tags);


};


YAFEL_NAMESPACE_CLOSE

#endif
