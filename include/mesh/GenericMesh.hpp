#ifndef __YAFEL_GENERIC_MESH_HPP
#define __YAFEL_GENERIC_MESH_HPP

/*
 * GenericMesh: Base class utilizing CRTP to define the interface for finite element
 * mesh data structures. It cannot be instantiated directly, but enables static
 * polymorphism for all functions taking a reference to a mesh.
 */

#include "yafel_globals.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include "mesh/MeshIterator.hpp"


#include <vector>

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

  // iterator begin() and end(). This should not be implemented in children of this class.
  MeshIterator<T> begin() const {
    return MeshIterator<T>(*this, 0);
  }
  MeshIterator<T> end() const {
    return MeshIterator<T>(*this, n_elements());
  }
  
  
  // Stuff for CRTP to work
  operator T&(){return static_cast<T&>(*this);}
  operator T const&() const {return static_cast<T const&>(*this);}
};

YAFEL_NAMESPACE_CLOSE

#endif
