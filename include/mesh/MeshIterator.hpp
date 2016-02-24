#ifndef __YAFEL_MESHITERATOR_HPP
#define __YAFEL_MESHITERATOR_HPP

#include "yafel_globals.hpp"
#include <exception>

/*
 * STL-like iterator to iterate over the elements in a mesh.
 * Maybe should make some distinction between element and node iterators,
 * forward and reverse iterators (though I can't think of a reason for both)?
 */

YAFEL_NAMESPACE_OPEN

template<typename MESHTYPE>
class MeshIterator {
public:
  using size_type = typename MESHTYPE::size_type;
  using coordinate_type = typename MESHTYPE::coordinate_type;
  using element_container = typename MESHTYPE::element_container;

  // Constructors
  MeshIterator(const MESHTYPE &M, size_type p) :
    mesh(M), ptr(p)
  {}


  // prefix increment
  MeshIterator<MESHTYPE> & operator++() {
    ++ptr; 
    return *this;
  }

  // postfix increment
  MeshIterator<MESHTYPE> operator++(int) {
    MeshIterator<MESHTYPE> ret(mesh,ptr); 
    ++ptr;
    return ret;
  }

  element_container operator*() const {
    return mesh.element(ptr);
  }
  
  bool operator<(const MeshIterator<MESHTYPE> &rhs) const {
    if(&mesh != &(rhs.mesh)) {
      // iterators do not point to the same mesh. Not cool bro...
      throw(std::invalid_argument("MeshIterator objects point to different meshes"));
    }
    return ptr<rhs.ptr;
  }

  bool operator==(const MeshIterator<MESHTYPE> &rhs) const {
    if(&mesh != &(rhs.mesh)) {
      // iterators do not point to the same mesh. Not cool bro...
      throw(std::invalid_argument("MeshIterator objects point to different meshes"));
    }
    return ptr==rhs.ptr;
  }
private:
  const MESHTYPE &mesh;
  size_type ptr;
  
};

YAFEL_NAMESPACE_CLOSE

#endif
