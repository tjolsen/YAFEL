#ifndef __YAFEL_GMSHMESH_HPP
#define __YAFEL_GMSHMESH_HPP


/* 
 * GmshMesh:
 * 
 * This data structure represents a mesh read from a gmsh .msh file.
 * Since these meshes have heterogeneous element types, it must explicitly
 * store the node locations, element connectivity, element types, and "tags"
 * (see gmsh documentation).
 * 
 */

#include "yafel_globals.hpp"
#include "mesh/GenericMesh.hpp"
#include <string>
#include <exception>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class GmshMesh : public GenericMesh<GmshMesh<NSD>, NSD> {

public:
  using coordinate_type = typename GenericMesh<GmshMesh<NSD>,NSD>::coordinate_type;
  using size_type = typename GenericMesh<GmshMesh<NSD>,NSD>::size_type;
  using element_container = typename GenericMesh<GmshMesh<NSD>,NSD>::element_container;

  std::vector<coordinate_type> _nodes;
  std::vector<element_container> _elements;
  std::vector<size_type> _element_type;
  std::vector<std::vector<size_type> > _element_tags;

  // Constructors
  GmshMesh(const std::string &fname) {
    
    
    
  }


  // Interface implementation
  inline size_type n_nodes() const {
    return _nodes.size();
  }
  
  // warning, this returns the number of ALL elements. Not just the useful ones.
  // Be smart when looping over the list of elements.
  inline size_type n_elements() const {
    return _elements.size();
  }

  inline coordinate_type node(size_type nodenum) const {
    return _nodes[nodenum];
  }

  inline element_container element(size_type elnum) const {
    return _elements[elnum];
  }
};

YAFEL_NAMESPACE_CLOSE

#endif
