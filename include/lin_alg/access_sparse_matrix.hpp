#ifndef _YAFEL_ACCESS_SPARSE_MATRIX
#define _YAFEL_ACCESS_SPARSE_MATRIX


#include "yafel_globals.hpp"
#include "lin_alg/sparse_matrix.hpp"


YAFEL_NAMESPACE_OPEN

template<typename T, typename dataType=double>
class access_sparse_matrix : public sparse_matrix<access_sparse_matrix<T,dataType>, dataType> {

public:
  using size_type = typename sparse_matrix<access_sparse_matrix<T,dataType>,dataType>::size_type;
  using value_type = typename sparse_matrix<access_sparse_matrix<T,dataType>,dataType>::value_type;
  using reference = typename sparse_matrix<access_sparse_matrix<T,dataType>,dataType>::reference;
  using triplet = typename sparse_matrix<access_sparse_matrix<T,dataType>,dataType>::triplet;

  // interface inherited from sparse_matrix
  size_type rows() const {return static_cast<const T&>(*this).rows();}
  size_type cols() const {return static_cast<const T&>(*this).cols();}
  size_type nnz() {return static_cast<T&>(*this).nnz();}
  const std::vector<triplet> & get_triplets() {return static_cast<T&>(*this).get_triplets();}
  std::vector<triplet> copy_triplets() {return static_cast<T&>(*this).copy_triplets();}

  operator T&() {return static_cast<T&>(*this);}
  operator T const&() {return static_cast<T const&>(*this);}


  /*
   * Interface for access matrices. Provides interface to access values/references in the matrix
   */
  value_type operator()(size_type i, size_type j) const {return static_cast<T const&>(*this)(i,j);}
  reference operator()(size_type i, size_type j) {return static_cast<T&>(*this)(i,j);}

};

YAFEL_NAMESPACE_CLOSE


#endif
