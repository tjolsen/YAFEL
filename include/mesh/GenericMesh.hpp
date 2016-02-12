#ifndef __YAFEL_GENERIC_MESH_HPP
#define __YAFEL_GENERIC_MESH_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/Tensor.hpp"

#include <vector>

/*
 * GenericMesh: Base class utilizing CRTP to define the interface for finite element
 * mesh data structures. 
 */

YAFEL_NAMESPACE_OPEN

template<typename T, unsigned NSD>
class GenericMesh {

public:
  using coordinate_type = Tensor<NSD,1,double>;
  using size_type = std::size_t;
  using element_container = std::vector<size_type>;

  // Functions that all meshes must support
  size_type n_nodes() const {return static_cast<T const&>(*this).n_nodes();}
  size_type n_elements() const {return static_cast<T const&>(*this).n_elements();}
  
  coordinate_type node(size_type nodenum) const {return static_cast<T const&>(*this).node(nodenum);}
  element_container element(size_type elnum) const {return static_cast<T const&>(*this).element(elnum);}


  // Stuff for CRTP to work
  operator T&(){return static_cast<T&>(*this);}
  operator T const&() const {return static_cast<T const&>(*this);}
};

YAFEL_NAMESPACE_CLOSE

#endif
