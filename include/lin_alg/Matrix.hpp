#ifndef __YAFEL_MATRIX_HPP
#define __YAFEL_MATRIX_HPP

#include "yafel_globals.hpp"
#include "lin_alg/MatrixExpression.hpp"
#include <cstdio>
#include <cstdlib>

YAFEL_NAMESPACE_OPEN

template<typename dataType=double>
class Matrix : public MatrixExpression<Matrix<dataType>,dataType> {

public:
  using container_type = typename MatrixExpression<Matrix<dataType>,dataType>::container_type;
  using value_type = typename MatrixExpression<Matrix<dataType>,dataType>::value_type;
  using size_type = typename MatrixExpression<Matrix<dataType>,dataType>::size_type;
  using reference = typename MatrixExpression<Matrix<dataType>,dataType>::reference;

  reference operator()(size_type i, size_type j) {return _data[linear_index(i,j)];}
  value_type operator()(size_type i, size_type j) const {
#ifndef _OPTIMIZED
    if(i<0 || i>= _rows) {
      fprintf(stderr, "Matrix(i,j): i out of bounds\n"
	      "i=%lu, _rows=%lu\n", i, _rows);
      exit(1);
    }
    if(j<0 || j>= _cols) {
      fprintf(stderr, "Matrix(i,j): j out of bounds\n"
	      "j=%lu, _cols=%lu\n", j, _cols);
      exit(1);
    }
#endif
    return _data[linear_index(i,j)];
  }
  size_type rows() const {return _rows;}
  size_type cols() const {return _cols;}
  
  /*
   * Constructors
   */
  // Square matrix of dimension N with default dataType construction
  Matrix(size_type N) : _data(N*N), _rows(N), _cols(N) {}

  // Matrix of dimensions M x N with default dataType construction
  Matrix(size_type r, size_type c) : _data(r*c), _rows(r), _cols(c) {}

  // Matrix of dimensions M x N with value dataType construction
  Matrix(size_type r, size_type c, dataType val) : _data(r*c, val), _rows(r), _cols(c) {}

  // Copy constructor
  //Matrix(const Matrix<dataType> &u) : _data(u._data), _rows(u._rows), _cols(u._cols) {}

  // Construct from MatrixExpression<T,dataType>
  template<typename T>
  Matrix(const MatrixExpression<T,dataType> &u) : _data(u.rows()*u.cols(), 0), 
						  _rows(u.rows()), _cols(u.cols()) {
    /*
    _rows = u.rows();
    _cols = u.cols();
    _data.resize(_rows*_cols);
    */
    for(size_type i=0; i<_rows; ++i) {
      for(size_type j=0; j<_cols; ++j) {
	(*this)(i,j) = u(i,j);
      }
    }
  }

  /*
   * Matrix-specific update operators. Cannot be used with MatrixExpressions
   * as left-hand argument.
   */
  template<typename T>
  Matrix<dataType>& operator+=(const MatrixExpression<T,dataType> &rhs) {
#ifndef _OPTIMIZED
    assert(rows() == rhs.rows() && cols() == rhs.cols() && "Matrix::operator+= dimension mismatch");
#endif
    for(size_type r=0; r<_rows; ++r) {
      for(size_type c=0; c<_cols; ++c) {
        (*this)(r,c) += rhs(r,c);
      }
    }

    return *this;
  }

  template<typename T>
  Matrix<dataType>& operator-=(const MatrixExpression<T,dataType> &rhs) {
#ifndef _OPTIMIZED
    assert(rows() == rhs.rows() && cols() == rhs.cols() && "Matrix::operator-= dimension mismatch");
#endif
    for(size_type r=0; r<_rows; ++r) {
      for(size_type c=0; c<_cols; ++c) {
        (*this)(r,c) -= rhs(r,c);
      }
    }
    
    return *this;
  }
  
  Matrix<dataType>& operator*=(dataType rhs) {
    for(size_type i=0; i<_data.size(); ++i) {
      _data[i] *= rhs;
    }

    return *this;
  }


private:
  container_type _data;
  size_type _rows;
  size_type _cols;

  /*
   * Define row-major linear indexing here to eliminate any explicit
   * dependence on this later in the code. Can be easily changed to column-major
   * if desired.
   */
  inline size_type linear_index(size_type i, size_type j) const {
    return i*_cols + j;
  }
};




YAFEL_NAMESPACE_CLOSE


#endif
