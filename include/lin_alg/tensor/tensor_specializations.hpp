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



YAFEL_NAMESPACE_CLOSE


#endif
