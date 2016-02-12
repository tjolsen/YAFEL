#ifndef _YAFEL_MESH_HPP
#define _YAFEL_MESH_HPP

#include <vector>
#include "yafel_globals.hpp"
#include "lin_alg/tensor/Tensor.hpp"

/*
 *
 * This header describes the interface for the Mesh class, which holds
 * node locations, element connectivity, element faces (tbd), additional "tags" (gmsh parlance)
 * associated with elements.
 *
 * 
 * The element_types vector holds a list of integers corresponding to the GMSH element type
 * (see online documentation for details)
 *
 */

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class Mesh {
  
public:
  using coordinate_type = Tensor<NSD,1,double>;
  using container_type = std::vector<coordinate_type>;
  using size_type = typename container_type::size_type;

  std::vector< coordinate_type > nodal_coords;
  std::vector<std::vector<size_type> > elements;
  std::vector< std::vector<size_type> > el_tags;
  std::vector<size_type> element_type;
  
  // ctors
  Mesh(const std::vector<coordinate_type> &nodes, 
       const std::vector<std::vector<size_type> > &elems, 
       const std::vector<size_type> & eltype,
       const std::vector<std::vector<size_type> > & _tags);

  Mesh(const std::vector<coordinate_type> &nodes, 
       const std::vector<std::vector<size_type> > &elems,
       const std::vector<size_type> & eltype);
  Mesh();

  // mesh characteristic size
  void compute_min_length();

  // getter functions
  inline size_type get_n_nodes() const {return this->n_nodes;}
  inline size_type get_n_elems() const {return this->n_elems;}
  inline double get_minLength() const {return this->minLength;}


  // Node reordering methods
  void reorder_rcm();


private:
  size_type n_elems;
  size_type n_nodes;
  double minLength;
  
  
};

YAFEL_NAMESPACE_CLOSE

#endif // end #ifndef _YAFEL_MESH_HP
