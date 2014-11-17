#ifndef _YAFEL_MESH_HPP
#define _YAFEL_MESH_HPP

#include <vector>
#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"

/*

This header describes the interface for the Mesh class, which holds
node locations, element connectivity, element faces (tbd), additional "tags" (gmsh parlance)
associated with elements.

The element_types vector holds a list of integers corresponding to the GMSH element type
(see online documentation for details)

 */

YAFEL_NAMESPACE_OPEN

class Mesh {
  
private:
  unsigned n_elems;
  unsigned n_nodes;
  
  public:
  std::vector< Vector > nodal_coords;
  std::vector< std::vector<int> > elements;
  std::vector< std::vector<int> > el_tags;
  std::vector<int> element_type;

  // set mesh nodes and elements
  Mesh();
  Mesh(const std::vector<Vector> &nodes, 
       const std::vector<std::vector<int> > &elems,
       const std::vector<int> & eltype);

  // set mesh nodes, elems, and tags
  Mesh(const std::vector<Vector> &nodes, 
       const std::vector<std::vector<int> > &elems, 
       const std::vector<int> & eltype,
       const std::vector<std::vector<int> > & _tags);
  
  inline unsigned get_n_nodes() {return this->n_nodes;}
  inline unsigned get_n_elems() {return this->n_elems;}
  
  // copy ctor -- use default
  //Mesh(const Mesh & m_in);

  //destructor -- use default
  //~Mesh();
  
};

YAFEL_NAMESPACE_CLOSE

#endif // end #ifndef _YAFEL_MESH_HPP
