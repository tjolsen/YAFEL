#ifndef _YAFEL_SPARSE_MATRIX_HPP
#define _YAFEL_SPARSE_MATRIX_HPP

#include "yafel_globals.hpp"

#include <vector>
#include <tuple>

YAFEL_NAMESPACE_OPEN

template<typename T, typename dataType=double>
class sparse_matrix {
  
public:
  using container_type = std::vector<dataType>;
  using size_type =  typename container_type::size_type;
  using value_type = typename container_type::value_type;
  using reference = typename container_type::reference;
  using triplet = std::tuple<size_type, size_type, value_type>;
  
  
  size_type rows() const { return static_cast<const T&>(*this).rows(); }
  size_type cols() const { return static_cast<const T&>(*this).cols(); }
  size_type nnz() {return static_cast<T&>(*this).nnz();}

  // moving this to the access_sparse_matrix CRTP-intermediate class.
  // This removes operator() from the construction_sparse_matrix class, taking capability away from sparse_coo
  //value_type operator()(size_type i, size_type j) const {return static_cast<const T&>(*this)(i,j);}
  
  operator T&() {return static_cast<T&>(*this);}
  operator T const&() {return static_cast<T const&>(*this);}

  const std::vector<triplet> & get_triplets() {return static_cast<T&>(*this).get_triplets();}
};

YAFEL_NAMESPACE_CLOSE

#endif
