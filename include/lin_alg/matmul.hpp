#ifndef __YAFEL_MATMUL_HPP
#define __YAFEL_MATMUL_HPP


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
 *
 * This can be compiled in either serial or parallel mode, determined by the presence of the
 * macro _YAFEL_PARALLEL_MATMUL
 *
 * This macro can/should be passed in as a compiler parameter via -D_YAFEL_PARALLEL_MATMUL
 */

#include "yafel_globals.hpp"
#include "lin_alg/MatrixExpression.hpp"
#include "lin_alg/Matrix.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <thread>

// AVX intrinsics
#include <immintrin.h>

YAFEL_NAMESPACE_OPEN

/* 
 * Hard-coded parameters to control when to terminate recursion and move into a
 * (potentially) tightly-optimized block multiplication kernel.
 */ 
const std::size_t recursion_cutoff = 128;

// limit thread-spawning depth to floor( log2(NProcessors) )
// This minimizes threading overhead while eliminating the need for
// more complex thread scheduling
constexpr std::size_t thread_depth_limit(std::size_t N) {
  return (N/2 == 0) ? 0 : 1 + thread_depth_limit(N/2);
}

const std::size_t max_threads = thread_depth_limit(std::thread::hardware_concurrency());

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
		    std::size_t p
#ifdef _YAFEL_PARALLEL_MATMUL
		    ,
		    std::size_t thread_depth
#endif
);


// forward declaration of matmul kernel function
template<typename dataType>
inline void matmul_kernel(dataType Ablock[recursion_cutoff][recursion_cutoff], 
		   dataType BblockT[recursion_cutoff][recursion_cutoff], 
		   dataType Cblock[recursion_cutoff][recursion_cutoff],
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
#ifdef _YAFEL_PARALLEL_MATMUL
  divconq_matmul(retmat, A, B, 0,0, 0,0, m,n,p,0);
#else
  divconq_matmul(retmat, A, B, 0,0, 0,0, m,n,p);
#endif
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
		    std::size_t p
#ifdef _YAFEL_PARALLEL_MATMUL
		    ,
		    std::size_t thread_depth
#endif
)
{
  
  //fprintf(stdout, "IL:%lu JL:%lu IR:%lu JR:%lu m:%lu n:%lu p:%lu\n",
  //	  ileft, jleft, iright, jright, m, n, p);

  
  // Base Case: Kernel for computing a small block using standard A_{ik} B_{kj} multiplication
  if(m<=recursion_cutoff &&
     n<=recursion_cutoff &&
     p<=recursion_cutoff) {

    // copy blocks into local memory blocks. May want to move these into a higher-up place later
    // if allocating them seems slow
    dataType Ablock[recursion_cutoff][recursion_cutoff]__attribute__((aligned(32)));
    dataType BblockT[recursion_cutoff][recursion_cutoff]__attribute__((aligned(32)));
    dataType Cblock[recursion_cutoff][recursion_cutoff]__attribute__((aligned(32)));
    
    for(std::size_t i=0; i<m; ++i) {
      for(std::size_t j=0; j<n; ++j) {
	Ablock[i][j] = A(ileft+i,jleft+j);
      }
    }
    // store the B block in transposed form, to make sequential accesses adjacent in memory
    // should also allow for future vectorization of matmul kernel
    for(std::size_t i=0; i<n; ++i) {
      for(std::size_t j=0; j<p; ++j) {
	BblockT[j][i] = B(iright+i,jright+j);
      }
    }

    // call the matrix multiplication kernel
    // a highly optimized version exists for intel core 4xxx or 5xxx processors (requires FMA instruction support)
    matmul_kernel(Ablock,BblockT,Cblock,m,n,p);
    
    // do update to large array
    for(std::size_t i=0; i<m; ++i) {
      for(std::size_t j=0; j<p; ++j) {
	C(ileft+i,jright+j) += Cblock[i][j];
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
    
#ifdef _YAFEL_PARALLEL_MATMUL
    if(thread_depth < max_threads) {
      /*
      std::thread t1([&](){divconq_matmul(C,A,B, iL1, jleft, iright, jright, 
					  m1, n, p, thread_depth+1);});
      std::thread t2([&](){divconq_matmul(C,A,B, iL2, jleft, iright, jright, 
					  m2, n, p, thread_depth+1);});
      */
      std::thread t1(divconq_matmul<T1,T2,dataType>, std::ref(C),std::ref(A),std::ref(B), iL1, jleft, iright, jright, 
		     m1, n, p, thread_depth+1);
      std::thread t2(divconq_matmul<T1,T2,dataType>, std::ref(C),std::ref(A),std::ref(B), iL2, jleft, iright, jright, 
		     m2, n, p, thread_depth+1);

      t1.join();
      t2.join();
    } else {
      divconq_matmul(C,A,B, iL1, jleft, iright, jright, m1, n, p, thread_depth);
      divconq_matmul(C,A,B, iL2, jleft, iright, jright, m2, n, p, thread_depth);
    }
#else
    divconq_matmul(C,A,B, iL1, jleft, iright, jright, m1, n, p);
    divconq_matmul(C,A,B, iL2, jleft, iright, jright, m2, n, p);
#endif
    
  }
  else if( p >= m && p >= n) {
    // split the "p" dimension
    std::size_t jR1 = jright;
    std::size_t p1 = p/2;
    std::size_t jR2 = jright + p1;
    std::size_t p2 = p - p1;
   
#ifdef _YAFEL_PARALLEL_MATMUL
    if(thread_depth < max_threads) {
      /*
      std::thread t1([&](){divconq_matmul(C,A,B, ileft, jleft, iright, jR1, 
					  m, n, p1, thread_depth+1);});
      std::thread t2([&](){divconq_matmul(C,A,B, ileft, jleft, iright, jR2, 
					  m, n, p2, thread_depth+1);});
      */
      std::thread t1(divconq_matmul<T1,T2,dataType>, std::ref(C),std::ref(A),std::ref(B), ileft, jleft, iright, jR1,
		     m, n, p1, thread_depth+1);
      std::thread t2(divconq_matmul<T1,T2,dataType>, std::ref(C),std::ref(A),std::ref(B), ileft, jleft, iright, jR2, 
		     m, n, p2, thread_depth+1);
      t1.join();
      t2.join();
    } else {
      divconq_matmul(C,A,B, ileft, jleft, iright, jR1, m, n, p1, thread_depth);
      divconq_matmul(C,A,B, ileft, jleft, iright, jR2, m, n, p2, thread_depth);
    }
#else
    divconq_matmul(C,A,B, ileft, jleft, iright, jR1, m, n, p1);
    divconq_matmul(C,A,B, ileft, jleft, iright, jR2, m, n, p2);
#endif
  }
  else if( n >= m && n >= p) {
    // split the "n" dimension
    std::size_t n1 = n/2;
    std::size_t n2 = n - n/2;
    
    std::size_t jLiR_1 = jleft;
    std::size_t jLiR_2 = jleft + n1;
    
#ifdef _YAFEL_PARALLEL_MATMUL
    divconq_matmul(C,A,B, ileft, jLiR_1, jLiR_1, jright, m, n1, p, thread_depth);
    divconq_matmul(C,A,B, ileft, jLiR_2, jLiR_2, jright, m, n2, p, thread_depth);
#else
    divconq_matmul(C,A,B, ileft, jLiR_1, jLiR_1, jright, m, n1, p);
    divconq_matmul(C,A,B, ileft, jLiR_2, jLiR_2, jright, m, n2, p);
#endif
  }
  
} // end divconq_matmul


template<typename dataType>
inline void 
matmul_kernel(dataType Ablock[recursion_cutoff][recursion_cutoff], 
		   dataType BblockT[recursion_cutoff][recursion_cutoff], 
		   dataType Cblock[recursion_cutoff][recursion_cutoff],
		   std::size_t m,
		   std::size_t n,
		   std::size_t p) {
  for(std::size_t i=0; i<m; ++i) {
    for(std::size_t j=0; j<p; ++j) {
      Cblock[i][j] = 0;
      
      // Sachith's inner loop. Can be adapated in partial specialization for dataType=double 
      // to use AVX FMA instructions?
      dataType s0 = 0;
      dataType s1 = 0;
      dataType s2 = 0;
      dataType s3 = 0;
      std::size_t nn = n & (-4); // 4*floor(n/4)
      for(std::size_t k=0; k<nn; k += 4) {
	s0 += Ablock[i][k+0]*BblockT[j][k+0];
	s1 += Ablock[i][k+1]*BblockT[j][k+1];
	s2 += Ablock[i][k+2]*BblockT[j][k+2];
	s3 += Ablock[i][k+3]*BblockT[j][k+3];
      }
      Cblock[i][j] += (s0 + s1 + s2 + s3);
      for (std::size_t k = nn; k<n; ++k){
	Cblock[i][j] += Ablock[i][k]*BblockT[j][k];
      }
    }
  }
  
  
} // end matmul_kernel


// Conditionally compile if processor supports FMA instructions
#ifdef __FMA__

template<>
//__attribute__((gnu_inline))
inline void 
matmul_kernel(double Ablock[recursion_cutoff][recursion_cutoff], 
		   double BblockT[recursion_cutoff][recursion_cutoff], 
		   double Cblock[recursion_cutoff][recursion_cutoff],
		   std::size_t m,
		   std::size_t n,
		   std::size_t p) {
  for(std::size_t i=0; i<m; ++i) {
    
    std::size_t pp = p & (-4);
    for(std::size_t j=0; j<pp; j += 4) {
      Cblock[i][j+0] = 0;
      Cblock[i][j+1] = 0;
      Cblock[i][j+2] = 0;
      Cblock[i][j+3] = 0;

      // For accumulating
      __m256d ymm0, ymm1, ymm2, ymm3;

      // For holding Ablock[i][k:k+3]
      __m256d ymm4, ymm5, ymm6, ymm7;
      
      // For holding BblockT[j+dj][k:k+3]
      __m256d ymm8, ymm9, ymm10, ymm11;
      
      ymm0 = _mm256_setzero_pd();
      ymm1 = _mm256_setzero_pd();
      ymm2 = _mm256_setzero_pd();
      ymm3 = _mm256_setzero_pd();

      std::size_t nn = 16*(n/16);// & (-16); // 16*floor(n/16)
      for(std::size_t k=0; k<nn; k += 16) {
	ymm4 = _mm256_load_pd(&(Ablock[i][k+0]));
	ymm5 = _mm256_load_pd(&(Ablock[i][k+4]));
	ymm6 = _mm256_load_pd(&(Ablock[i][k+8]));
	ymm7 = _mm256_load_pd(&(Ablock[i][k+12]));

	//---------- dk=0 -----------
	ymm8 = _mm256_load_pd(&(BblockT[j+0][k+0]));
	ymm9 = _mm256_load_pd(&(BblockT[j+1][k+0]));
	ymm10 = _mm256_load_pd(&(BblockT[j+2][k+0]));
	ymm11 = _mm256_load_pd(&(BblockT[j+3][k+0]));

	ymm0 = _mm256_fmadd_pd(ymm4, ymm8, ymm0);
	ymm1 = _mm256_fmadd_pd(ymm4, ymm9, ymm1);
	ymm2 = _mm256_fmadd_pd(ymm4, ymm10, ymm2);
	ymm3 = _mm256_fmadd_pd(ymm4, ymm11, ymm3);

	//---------- dk=4 -----------
	ymm8 = _mm256_load_pd(&(BblockT[j+0][k+4]));
	ymm9 = _mm256_load_pd(&(BblockT[j+1][k+4]));
	ymm10 = _mm256_load_pd(&(BblockT[j+2][k+4]));
	ymm11 = _mm256_load_pd(&(BblockT[j+3][k+4]));

	ymm0 = _mm256_fmadd_pd(ymm5, ymm8, ymm0);
	ymm1 = _mm256_fmadd_pd(ymm5, ymm9, ymm1);
	ymm2 = _mm256_fmadd_pd(ymm5, ymm10, ymm2);
	ymm3 = _mm256_fmadd_pd(ymm5, ymm11, ymm3);

	
	//---------- dk=8 -----------
	ymm8 = _mm256_load_pd(&(BblockT[j+0][k+8]));
	ymm9 = _mm256_load_pd(&(BblockT[j+1][k+8]));
	ymm10 = _mm256_load_pd(&(BblockT[j+2][k+8]));
	ymm11 = _mm256_load_pd(&(BblockT[j+3][k+8]));

	ymm0 = _mm256_fmadd_pd(ymm6, ymm8, ymm0);
	ymm1 = _mm256_fmadd_pd(ymm6, ymm9, ymm1);
	ymm2 = _mm256_fmadd_pd(ymm6, ymm10, ymm2);
	ymm3 = _mm256_fmadd_pd(ymm6, ymm11, ymm3);

	//---------- dk=12 -----------
	ymm8 = _mm256_load_pd(&(BblockT[j+0][k+12]));
	ymm9 = _mm256_load_pd(&(BblockT[j+1][k+12]));
	ymm10 = _mm256_load_pd(&(BblockT[j+2][k+12]));
	ymm11 = _mm256_load_pd(&(BblockT[j+3][k+12]));

	ymm0 = _mm256_fmadd_pd(ymm7, ymm8, ymm0);
	ymm1 = _mm256_fmadd_pd(ymm7, ymm9, ymm1);
	ymm2 = _mm256_fmadd_pd(ymm7, ymm10, ymm2);
	ymm3 = _mm256_fmadd_pd(ymm7, ymm11, ymm3);


      }

      // horizontal add to accumulate dot prods
      Cblock[i][j+0] += ( ((double*)&ymm0)[0]+((double*)&ymm0)[1]+((double*)&ymm0)[2]+((double*)&ymm0)[3] );
      Cblock[i][j+1] += ( ((double*)&ymm1)[0]+((double*)&ymm1)[1]+((double*)&ymm1)[2]+((double*)&ymm1)[3] );
      Cblock[i][j+2] += ( ((double*)&ymm2)[0]+((double*)&ymm2)[1]+((double*)&ymm2)[2]+((double*)&ymm2)[3] );
      Cblock[i][j+3] += ( ((double*)&ymm3)[0]+((double*)&ymm3)[1]+((double*)&ymm3)[2]+((double*)&ymm3)[3] );

      for (std::size_t k = nn; k<n; ++k){
	Cblock[i][j+0] += Ablock[i][k]*BblockT[j+0][k];
	Cblock[i][j+1] += Ablock[i][k]*BblockT[j+1][k];
	Cblock[i][j+2] += Ablock[i][k]*BblockT[j+2][k];
	Cblock[i][j+3] += Ablock[i][k]*BblockT[j+3][k];
      }
    } // end j<pp loop
    for (std::size_t j = pp; j<p; ++j) {
      Cblock[i][j] = 0;

      double s0 = 0;
      double s1 = 0;
      double s2 = 0;
      double s3 = 0;
      std::size_t nn = n &(-4);
      for(std::size_t k=0; k<nn; k+= 4) {
	s0 += Ablock[i][k+0]*BblockT[j][k+0];
	s1 += Ablock[i][k+1]*BblockT[j][k+1];
	s2 += Ablock[i][k+2]*BblockT[j][k+2];
	s3 += Ablock[i][k+3]*BblockT[j][k+3];
      }
      Cblock[i][j] += (s0 + s1 + s2 + s3);
      for(std::size_t k=nn; k<n; ++k) {
	Cblock[i][j] += Ablock[i][k]*BblockT[j][k];
      }
    }
  }
  
  
} // end matmul_kernel

#endif //end #ifdef __FMA__

YAFEL_NAMESPACE_CLOSE

#endif
