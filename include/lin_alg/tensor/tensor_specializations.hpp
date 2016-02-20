#ifndef __YAFEL_TENSOR_SPECIALIZATIONS_HPP
#define __YAFEL_TENSOR_SPECIALIZATIONS_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"


YAFEL_NAMESPACE_OPEN


// specialize operator* for rank-2 * rank-1 tensor expressions (contract 1)
template<typename T1, typename T2, unsigned DIM, typename dataType>
TensorContraction<T1,T2,DIM,2,1,1,dataType>
operator*(const TensorExpression<T1,DIM,2,dataType> &lhs,
	  const TensorExpression<T2,DIM,1,dataType> &rhs) {
  //return TensorContraction<T1,T2,DIM,2,1,1,dataType>(lhs,rhs);
  return contract<1>(lhs,rhs);
}

// specialize operator* for rank-2 * rank-2 tensor expressions (contract 1)
template<typename T1, typename T2, unsigned DIM, typename dataType>
TensorContraction<T1,T2,DIM,2,2,1,dataType>
operator*(const TensorExpression<T1,DIM,2,dataType> &lhs,
	  const TensorExpression<T2,DIM,2,dataType> &rhs) {
  //return TensorContraction<T1,T2,DIM,2,2,1,dataType>(lhs,rhs);
  return contract<1>(lhs,rhs);
}

// specialize operator* for rank-4 * rank-2 tensor expressions (contract 2)
template<typename T1, typename T2, unsigned DIM, typename dataType>
TensorContraction<T1,T2,DIM,4,2,2,dataType>
operator*(const TensorExpression<T1,DIM,4,dataType> &lhs,
	  const TensorExpression<T2,DIM,2,dataType> &rhs) {
  //return TensorContraction<T1,T2,DIM,4,2,2,dataType>(lhs,rhs);
  return contract<2>(lhs,rhs);
}


// determinant of rank=2, dim=2,3 tensors
template<typename T, typename dataType>
dataType det(const TensorExpression<T,2,2,dataType> &A) {
  return A(0,0)*A(1,1) - A(1,0)*A(0,1);
}
template<typename T, typename dataType>
dataType det(const TensorExpression<T,3,2,dataType> &A) {
  return A(0,0)*(A(1,1)*A(2,2) - A(2,1)*A(1,2))
    - A(0,1)*(A(1,0)*A(2,2) - A(2,0)*A(1,2))
    + A(0,2)*(A(1,0)*A(2,1) - A(2,0)*A(1,1));
}


// inverse of rank=2, dim=2,3 tensors
template<typename T, typename dataType>
Tensor<2,2,dataType> inv(const TensorExpression<T,2,2,dataType> &A) {
  
  dataType detA = det(A);
  Tensor<2,2,dataType> Ainv;
  
  Ainv(0,0) = A(1,1)/detA;
  Ainv(0,1) = -A(0,1)/detA;
  Ainv(1,0) = -A(1,0)/detA;
  Ainv(1,1) = A(0,0)/detA;
  
  return Ainv;
}

template<typename T, typename dataType>
Tensor<3,2,dataType> inv(const TensorExpression<T,3,2,dataType> &A) {
  
  dataType detA = det(A);
  Tensor<3,2,dataType> Ainv;
  
  Ainv(0,0) =  (A(1,1)*A(2,2) - A(2,1)*A(1,2))/detA;
  Ainv(0,1) = -(A(0,1)*A(2,2) - A(0,2)*A(2,1))/detA;
  Ainv(0,2) =  (A(0,1)*A(1,2) - A(0,2)*A(1,1))/detA;
  Ainv(1,0) = -(A(1,0)*A(2,2) - A(1,2)*A(2,0))/detA;
  Ainv(1,1) =  (A(0,0)*A(2,2) - A(0,2)*A(2,0))/detA;
  Ainv(1,2) = -(A(0,0)*A(1,2) - A(0,2)*A(1,0))/detA;
  Ainv(2,0) =  (A(1,0)*A(2,1) - A(1,1)*A(2,0))/detA;
  Ainv(2,1) = -(A(0,0)*A(2,1) - A(0,1)*A(2,0))/detA;
  Ainv(2,2) =  (A(0,0)*A(1,1) - A(0,1)*A(1,0))/detA;
  return Ainv;
}

YAFEL_NAMESPACE_CLOSE


#endif
