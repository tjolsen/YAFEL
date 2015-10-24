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

#include <cassert>
#include <cstdlib>
#include <iostream>

YAFEL_NAMESPACE_OPEN

/* 
 * Hard-coded parameters to control when to terminate recursion and move into a
 * (potentially) tightly-optimized block multiplication kernel.
 */ 
const std::size_t recursion_cutoff = 32;




//forward function declaration
template<typename T1, typename T2, typename dataType>
void divconq_matmul(Matrix<dataType> & C, 
		    const MatrixExpression<T1, dataType> & A,
		    const MatrixExpression<T2, dataType> & B,
		    std::size_t ileft,
		    std::size_t jleft,
		    std::size_t iright,
		    std::size_t jright,
		    std::size_t m,
		    std::size_t n,
		    std::size_t p);


// top-level matrix multiplication function
template<typename T1, typename T2, typename dataType>
Matrix<dataType> matmul(const MatrixExpression<T1,dataType> &A,
			const MatrixExpression<T2,dataType> &B) {

#ifndef _OPTIMIZED
  // Ensure proper dimensions
  assert(A.cols() == B.rows() && "Error: Matmul: MatrixExpression dimension mismatch.");
#endif

  std::size_t m,n,p;
  m = A.rows();
  n = B.rows();
  p = B.cols();

  // allocate return matrix
  Matrix<dataType> retmat(m, p, 0);

  // call the recursive matrix multiply
  divconq_matmul(retmat, A, B, 0,0, 0,0, m,n,p);

  return retmat;
}


/*
 * Divide-and-conquer recursive subroutine to compute C = A*B
 * {ileft,jleft} = upper-left corner of block of A being considered
 * {iright,jright} = upper-left corner of block of B being considered
 * m = number of rows in A
 * n = number of cols in A AND number of rows in B
 * p = number of cols in B
 */
template<typename T1, typename T2, typename dataType>
void divconq_matmul(Matrix<dataType> & C, 
		    const MatrixExpression<T1, dataType> & A,
		    const MatrixExpression<T2, dataType> & B,
		    std::size_t ileft,
		    std::size_t jleft,
		    std::size_t iright,
		    std::size_t jright,
		    std::size_t m,
		    std::size_t n,
		    std::size_t p)
{
  
  //fprintf(stdout, "IL:%lu JL:%lu IR:%lu JR:%lu m:%lu n:%lu p:%lu\n",
  //	 ileft, jleft, iright, jright, m, n, p);
  
  
  // Base Case: Kernel for computing a small block using standard A_{ik} B_{kj} multiplication
  if(m<=recursion_cutoff &&
     n<=recursion_cutoff &&
     p<=recursion_cutoff) {

    // copy blocks into local memory blocks. May want to move these into a higher-up place later
    // if allocating them seems slow
    dataType Ablock[recursion_cutoff][recursion_cutoff];
    dataType BblockT[recursion_cutoff][recursion_cutoff];
    for(std::size_t i=0; i<m; ++i) {
      for(std::size_t j=0; j<n; ++j) {
	Ablock[i][j] = A(i,j);
      }
    }
    // store the B block in transposed form, to make sequential accesses adjacent in memory
    // should also allow for future vectorization of matmul kernel
    for(std::size_t j=0; j<p; ++j) {
      for(std::size_t i=0; i<n; ++i) {
	BblockT[j][i] = B(i,j);
      }
    }

    for(std::size_t i=0; i<m; ++i) {
      for(std::size_t j=0; j<p; ++j) {
	for(std::size_t k=0; k<n; ++k) {
	  
	  C(ileft+i, jright+j) += Ablock[i][k]*BblockT[j][k];
	  
	}
      }
    }
    
    return;
  }
  

  // Recursion: Note, these cases ARE mutually exclusive. One of {m,n,p} must be the largest of the three
  if( m >= n && m >= p ) {
    //split the "m" dimension
    
    std::size_t iL1 = ileft;
    std::size_t m1 = m/2;
    std::size_t iL2 = ileft + m1;
    std::size_t m2 = m-m1;
    
    divconq_matmul(C,A,B, iL1, jleft, iright, jright, m1, n, p);
    divconq_matmul(C,A,B, iL2, jleft, iright, jright, m2, n, p);
  }
  else if( n >= m && n >= p) {
    // split the "n" dimension
    std::size_t n1 = n/2;
    std::size_t n2 = n - n/2;
    
    std::size_t jLiR_1 = jleft;
    std::size_t jLiR_2 = jleft + n1;
    
    divconq_matmul(C,A,B, ileft, jLiR_1, jLiR_1, jright, m, n1, p);
    divconq_matmul(C,A,B, ileft, jLiR_2, jLiR_2, jright, m, n2, p);
  }
  else if( p >= m && p >= n) {
    // split the "p" dimension
    
    std::size_t jR1 = jright;
    std::size_t p1 = p/2;
    std::size_t jR2 = jright + p1;
    std::size_t p2 = p - p1;
    
    divconq_matmul(C,A,B, ileft, jleft, iright, jR1, m, n, p1);
    divconq_matmul(C,A,B, ileft, jleft, iright, jR2, m, n, p2);
  }

}


YAFEL_NAMESPACE_CLOSE

#endif
