#ifndef _YAFEL_CSR_SPMV
#define _YAFEL_CSR_SPMV


/*
 * Sparse matrix-vector multiplication subroutines for use with matrices
 * in compressed sparse row format. This will only be implemented for Vector<dataType>
 * objects rather than VectorExpressions (implicit conversion via ctors will work)
 * due to the lack of data reuse in the spmv algorithm (unlike, for example, dense matmul).
 *
 * The algorithm returns a Vector<dataType> object.
 */

#include "yafel_globals.hpp"
#include "Vector.hpp"
#include "sparse_csr.hpp"
#include "csr_sparsity_pattern.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>

YAFEL_NAMESPACE_OPEN

template<typename dataType>
Vector<dataType> csr_spmv(const sparse_csr<dataType> &A, 
			  const Vector<dataType> &x) 
{

#ifndef _OPTIMIZED
  assert(A.cols() == x.size() && "Error: csr_spmv dimension mismatch");
#endif
  
  Vector<dataType> b(A.rows(), 0);

#ifdef _OPENMP 
#pragma omp parallel for schedule(static,128)
#endif
  for(std::size_t row = 0; row<A.rows(); ++row) {

    std::size_t idxmin, idxmax;
    idxmin = A.row_ptr[row];
    idxmax = A.row_ptr[row+1];
    
    //accumulator
    dataType accum = 0;
    for(std::size_t idx=idxmin; idx<idxmax; ++idx) {
      std::size_t col = A.col_index[idx];
      accum += A._data[idx]*x(col);
    }
    b(row) = accum;
  }
  
  
  return b;
}



YAFEL_NAMESPACE_CLOSE


#endif
