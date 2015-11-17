#ifndef __YAFEL_CONSTRUCTION_SPARSE_MATRIX
#define __YAFEL_CONSTRUCTION_SPARSE_MATRIX

#include "yafel_globals.hpp"
#include "lin_alg/sparse_matrix.hpp"

YAFEL_NAMESPACE_OPEN

template<typename T, typename dataType=double>
class construction_sparse_matrix : public sparse_matrix<construction_sparse_matrix<T,dataType>, dataType> {

public:
  using size_type = typename sparse_matrix<construction_sparse_matrix<T,dataType>,dataType>::size_type;
  using value_type = typename sparse_matrix<construction_sparse_matrix<T,dataType>,dataType>::value_type;
  using reference = typename sparse_matrix<construction_sparse_matrix<T,dataType>,dataType>::reference;
  using triplet = typename sparse_matrix<construction_sparse_matrix<T,dataType>,dataType>::triplet;
  

  //interface inherited from sparse_matrix (must be implemented here like this to get static polymorphism right)
  size_type rows() const { return static_cast<const T&>(*this).rows(); }
  size_type cols() const { return static_cast<const T&>(*this).cols(); }
  size_type nnz() {return static_cast<T&>(*this).nnz();}
  //value_type operator()(size_type i, size_type j) const {return static_cast<const T&>(*this)(i,j);}
  const std::vector<triplet> & get_triplets() {return static_cast<T&>(*this).get_triplets();}
  
  operator T&() {return static_cast<T&>(*this);}
  operator T const&() {return static_cast<T const&>(*this);}


  /*
   * new interface for construction_sparse_matrix
   * This allows children of this class to efficiently implement a method to add elements
   * incrementally to construct a sparse matrix
   */
  void add(size_type i, size_type j, value_type val) {static_cast<T&>(*this).add(i,j,val);}

  // Suggest a new size for the underlying data storage to speed up incremental assembly
  // Will not do anything if it would result in the destruction of data (ie, N < data.size())
  void preallocate(size_type N) {static_cast<T&>(*this).preallocate(N);}
};

YAFEL_NAMESPACE_CLOSE

#endif
