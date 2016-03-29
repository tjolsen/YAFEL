#ifndef _YAFEL_MVMUL_HPP
#define _YAFEL_MVMUL_HPP

/*
 * Dense MatrixExpression-VectorExpression multiplication.
 * A Vector<dataType> object will be returned from this subroutine
 */
#include "yafel_globals.hpp"
#include "lin_alg/MatrixExpression.hpp"
#include "lin_alg/VectorExpression.hpp"
#include "lin_alg/Vector.hpp"

#include <exception>
#include <cstdlib>
#include <iostream>

#include <thread>

//AVX Intrinsics
#include <immintrin.h>

YAFEL_NAMESPACE_OPEN

/*
 * Naive implementation. Just get something working immediately.
 *
 * Should implement things like matmul later, but this is likely a memory-bandwidth-constrained problem
 */
template<typename MT, typename VT, typename dataType>
Vector<dataType> mvmul(const MatrixExpression<MT,dataType> &A,
                       const VectorExpression<VT,dataType> &x) {

#ifndef _OPTIMIZED
  //ensure proper dimensions
  if(A.cols() != x.size()) {
    throw std::length_error("mvmul dimension mismatch");
  }
#endif
  
  Vector<dataType> b(A.rows());
  
  for(std::size_t i=0; i<A.rows(); ++i) {
    dataType sum(0);
    for(std::size_t j=0; j<A.cols(); ++j) {
      sum += A(i,j)*x(j);
    }
    b(i) = sum;
  }
  
  return b;
}


YAFEL_NAMESPACE_CLOSE

#endif
