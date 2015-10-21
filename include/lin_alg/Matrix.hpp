#ifndef __YAFEL_MATRIX_HPP
#define __YAFEL_MATRIX_HPP

#include "yafel_globals.hpp"
#include "lin_alg/MatrixExpression.hpp"

YAFEL_NAMESPACE_OPEN

template<typename dataType=double>
class Matrix : public MatrixExpression<Matrix<dataType>,dataType> {

private:
  typedef typename MatrixExpression<Matrix<dataType>,dataType>::container_type container_type;
  typedef typename MatrixExpression<Matrix<dataType>,dataType>::value_type value_type;
  typedef typename MatrixExpression<Matrix<dataType>,dataType>::size_type size_type;
  typedef typename MatrixExpression<Matrix<dataType>,dataType>::reference reference;
  container_type _data;
  size_type _rows;
  size_type _cols;
  
  /* 
   * Define row-major linear indexing.
   * Auxiliary function is used to remove any hidden explicit dependence on row/col-major
   * ordering elsewhere in code. Can be switched easily to column-major by 
   */
  inline size_type linear_index(i,j) const {return i*cols + j;}

public:
  reference operator()(size_type i, size_type j) {return _data[linear_index(i,j)];}
  value_type operator()(size_type i, size_type j) const {return _data[linear_index(i,j)];}

  /*
   * Constructors
   */
  // Square matrix of dimension N with default dataType construction
  Matrix(size_type N) : _data(N), _rows(N), _cols(N) {}

  // Matrix of dimensions M x N with default dataType construction
  Matrix(size_type r, size_type c) : _data(r*c), _rows(r), _cols(c) {}

  // Construct from MatrixExpression<T,dataType>
  template<typename T>
  Matrix(const MatrixExpression<T,dataType> &u) {
    
  }
};




YAFEL_NAMESPACE_CLOSE


#endif
