#ifndef _YAFEL_MATMUL_HPP
#define _YAFEL_MATMUL_HPP


/*
 * Matrix-matrix multiplication subroutines for use with MatrixExpressions.
 * A full Matrix<dataType> object is to be returned from the subroutines from
 * the product of two MatrixExpression<T1,dataType> and MatrixExpression<T2,dataType> objects.
 *
 *
 * The algorithm will use a cache-oblivious divide-and-conquer matrix multiplication algorithm
 * outlined on the cache-oblivious matrix multiplication wikipedia page. The original
 * algorithm is formalized in Harald Prokop's 1999 MIT SM thesis.
 * https://en.wikipedia.org/wiki/Cache-oblivious_matrix_multiplication
 * http://supertech.csail.mit.edu/papers/Prokop99.pdf
 */

#include "yafel_globals.hpp"
#include "lin_alg/MatrixExpression.hpp"
#include "lin_alg/Matrix.hpp"

YAFEL_NAMESPACE_OPEN








YAFEL_NAMESPACE_CLOSE

#endif
