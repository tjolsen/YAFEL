#ifndef __YAFEL_OPERATORS_HPP
#define __YAFEL_OPERATORS_HPP

#include "yafel_globals.hpp"

#include "lin_alg/Vector.hpp"
#include "lin_alg/VectorExpression.hpp"

#include "lin_alg/Matrix.hpp"
#include "lin_alg/MatrixExpression.hpp"
#include "lin_alg/matmul.hpp"
#include "lin_alg/mvmul.hpp"

#include "lin_alg/access_sparse_matrix.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/sparse_bcsr.hpp"
#include "lin_alg/csr_spmv.hpp"
#include "lin_alg/bcsr_spmv.hpp"



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
template<typename T1, typename T2, typename dataType, typename = typename std::enable_if<std::is_fundamental<T2>::value,T2>::type>
VectorScaled<T1,dataType> operator*(const VectorExpression<T1,dataType> &v, T2 a) {
  return VectorScaled<T1,dataType>(v, a);
}
template<typename T1, typename T2, typename dataType, typename = typename std::enable_if<std::is_fundamental<T2>::value,T2>::type>
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
template<typename T1, typename T2, typename dataType, typename = typename std::enable_if<std::is_fundamental<T2>::value,T2>::type>
MatrixScaled<T1,dataType> operator*(const MatrixExpression<T1,dataType> &v, T2 a) {
  return MatrixScaled<T1,dataType>(v,a);
}
template<typename T1, typename T2, typename dataType, typename = typename std::enable_if<std::is_fundamental<T2>::value,T2>::type>
MatrixScaled<T1,dataType> operator*(T2 a, const MatrixExpression<T1,dataType> &v) {
  return MatrixScaled<T1,dataType>(v,a);
}
//----------------------------------------------------------------------------------------
template<typename MT, typename VT, typename dataType>
Vector<dataType> operator*(const MatrixExpression<MT,dataType> &lhs,
                           const VectorExpression<VT,dataType> &rhs) {
  return mvmul(lhs, rhs);
}


//----------------------------------------------------------------------------------------
/*
 * Sparse matrix operators
 */
//----------------------------------------------------------------------------------------
template<typename T, typename dataType>
Vector<dataType> operator*(const access_sparse_matrix<T,dataType> &A,
                           const Vector<dataType> &x) {
  return static_cast<T const&>(A)*x;
}

template<unsigned BLOCK, typename dataType>
Vector<dataType> operator*(const sparse_bcsr<BLOCK, dataType> &A, const Vector<dataType> &x) {
  return bcsr_spmv(A,x);
}

template<typename dataType>
Vector<dataType> operator*(const sparse_csr<dataType> &A, const Vector<dataType> &x) {
  return csr_spmv(A,x);
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
