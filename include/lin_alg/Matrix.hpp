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
  inline size_type linear_index(size_type i, size_type j) const {return i*_cols + j;}

public:
  reference operator()(size_type i, size_type j) {return _data[linear_index(i,j)];}
  value_type operator()(size_type i, size_type j) const {return _data[linear_index(i,j)];}
  size_type rows() const {return _rows;}
  size_type cols() const {return _cols;}
  
  /*
   * Constructors
   */
  // Square matrix of dimension N with default dataType construction
  Matrix(size_type N) : _data(N), _rows(N), _cols(N) {}

  // Matrix of dimensions M x N with default dataType construction
  Matrix(size_type r, size_type c) : _data(r*c), _rows(r), _cols(c) {}

  // Matrix of dimensions M x N with value dataType construction
  Matrix(size_type r, size_type c, dataType val) : _data(r*c, val), _rows(r), _cols(c) {}

  // Construct from MatrixExpression<T,dataType>
  template<typename T>
  Matrix(const MatrixExpression<T,dataType> &u) {
    _rows = u.rows();
    _cols = u.cols();
    _data.resize(_rows*_cols);
    for(size_type i=0; i<_rows; ++i) {
      for(size_type j=0; j<_cols; ++j) {
	(*this)(i,j) = u(i,j);
      }
    }
  }
};




YAFEL_NAMESPACE_CLOSE


#endif
