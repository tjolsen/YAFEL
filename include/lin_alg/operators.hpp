#ifndef __YAFEL_OPERATORS_HPP
#define __YAFEL_OPERATORS_HPP

#include "yafel_globals.hpp"

#include "lin_alg/Vector.hpp"
#include "lin_alg/VectorExpression.hpp"

#include "lin_alg/Matrix.hpp"
#include "lin_alg/MatrixExpression.hpp"

#include "lin_alg/matmul.hpp"

YAFEL_NAMESPACE_OPEN

//-------------------------------------------------------------------------
/*
 * Operator overloading: Define operators (+,-,*,...) for VectorExpressions
 */
//-------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
VectorSum<T1,T2,dataType> operator+(const VectorExpression<T1,dataType> &lhs,
				    const VectorExpression<T2,dataType> &rhs) {
  return VectorSum<T1,T2,dataType>(lhs,rhs);
}

//-------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
VectorDifference<T1,T2,dataType> operator-(const VectorExpression<T1,dataType> &lhs,
				    const VectorExpression<T2,dataType> &rhs) {
  return VectorDifference<T1,T2,dataType>(lhs,rhs);
}

//-------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
VectorScaled<T1,dataType> operator*(const VectorExpression<T1,dataType> &v, T2 a) {
  return VectorScaled<T1,dataType>(v, a);
}
template<typename T1, typename T2, typename dataType>
VectorScaled<T1,dataType> operator*(T2 a, const VectorExpression<T1,dataType> &v) {
  return VectorScaled<T1,dataType>(v, a);
}

//-------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
bool operator==(const VectorExpression<T1,dataType> &lhs,
		const VectorExpression<T2,dataType> &rhs) {
  if(lhs.size() != rhs.size()) {
    return false;
  }
  
  for(typename VectorExpression<T1,dataType>::size_type i=0; i<lhs.size(); ++i) {
    if(lhs(i) != rhs(i)) {
      return false;
    }
  }
  
  return true;
}

//----------------------------------------------------------------------------------------
/*
 * Operator overloading: Define operators (+, -, *, ...) for MatrixExpressions
 */
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
MatrixSum<T1,T2,dataType> operator+(const MatrixExpression<T1,dataType> &lhs,
				    const MatrixExpression<T2,dataType> &rhs) {
  return MatrixSum<T1,T2,dataType>(lhs,rhs);
}
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
MatrixDifference<T1,T2,dataType> operator-(const MatrixExpression<T1,dataType> &lhs,
					   const MatrixExpression<T2,dataType> &rhs) {
  return MatrixDifference<T1,T2,dataType>(lhs,rhs);
}
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
Matrix<dataType> operator*(const MatrixExpression<T1,dataType> &lhs,
			   const MatrixExpression<T2,dataType> &rhs) {
  return matmul(lhs,rhs);
}
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
MatrixScaled<T1,dataType> operator*(const MatrixExpression<T1,dataType> &v, T2 a) {
  return MatrixScaled<T1,dataType>(v,a);
}
template<typename T1, typename T2, typename dataType>
MatrixScaled<T1,dataType> operator*(T2 a, const MatrixExpression<T1,dataType> &v) {
  return MatrixScaled<T1,dataType>(v,a);
}

//----------------------------------------------------------------------------------------
/*
 * Some relational operators that may be userful
 */
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
bool operator==(const MatrixExpression<T1,dataType> &lhs,
		const MatrixExpression<T2,dataType> &rhs) {
  
  if(lhs.rows()!=rhs.rows() || lhs.cols()!=rhs.cols()) {
    return false;
  }

  for(typename MatrixExpression<T1,dataType>::size_type i=0; i<lhs.rows(); ++i) {
    for(typename MatrixExpression<T1,dataType>::size_type j=0; j<lhs.cols(); ++j) {
      if(lhs(i,j) != rhs(i,j)) {
	return false;
      }
    }
  }
  
  return true;
}


YAFEL_NAMESPACE_CLOSE

#endif
