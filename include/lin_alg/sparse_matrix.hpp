#ifndef _YAFEL_SPARSE_MATRIX_HPP
#define _YAFEL_SPARSE_MATRIX_HPP

#include "yafel_globals.hpp"

#include <vector>

YAFEL_NAMESPACE_OPEN

template<typename T, typename dataType>
class sparse_matrix {
  
public:
  typedef std::vector<dataType> container_type;
  typedef typename container_type::size_type size_type;
  typedef typename container_type::value_type value_type;
  typedef typename container_type::reference reference;
  
  size_type rows() const { return static_cast<const T&>(*this).rows(); }
  size_type cols() const { return static_cast<const T&>(*this).cols(); }
  value_type operator()(size_type i, size_type j) const {return static_cast<const T&>(*this)(i,j);}
  
  operator T&() {return static_cast<T&>(*this);}
  operator T const&() {return static_cast<T const&>(*this);}
};

YAFEL_NAMESPACE_CLOSE

#endif
