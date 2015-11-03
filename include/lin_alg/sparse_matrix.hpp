#ifndef _YAFEL_SPARSE_MATRIX_HPP
#define _YAFEL_SPARSE_MATRIX_HPP

#include "yafel_globals.hpp"
#include "lin_alg/sparse_utils.hpp"

#include <vector>
#include <tuple>

YAFEL_NAMESPACE_OPEN

template<typename T, typename dataType=double>
class sparse_matrix {
  
public:
  /*
  typedef std::vector<dataType> container_type;
  typedef typename container_type::size_type size_type;
  typedef typename container_type::value_type value_type;
  typedef typename container_type::reference reference;
  */
  
  typedef std::size_t size_type;
  typedef dataType value_type;
  typedef dataType& reference;

  
  
  size_type rows() const { return static_cast<const T&>(*this).rows(); }
  size_type cols() const { return static_cast<const T&>(*this).cols(); }
  size_type nnz() {return static_cast<T&>(*this).nnz();}
  value_type operator()(size_type i, size_type j) const {return static_cast<const T&>(*this)(i,j);}
  
  operator T&() {return static_cast<T&>(*this);}
  operator T const&() {return static_cast<T const&>(*this);}

  std::vector<triplet> get_triplets() return {return static_cast<T&>(*this).get_triplets();}
};

YAFEL_NAMESPACE_CLOSE

#endif
