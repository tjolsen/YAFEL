#ifndef __YAFEL_CG_SOLVER_HPP
#define __YAFEL_CG_SOLVER_HPP

/*
 * Implementation of the conjugate-gradient iterative solver
 * for lienar systems with positive-definite, hermitian matrices.
 * The implementation follows the algorithm given in the wikipedia
 * article, which I believe follows Trefethen.
 */ 

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/access_sparse_matrix.hpp"
#include "lin_alg/operators.hpp"
#include <iostream>

YAFEL_NAMESPACE_OPEN

// not sure if I should use double or dataType for this.
// practically speaking, cg only makes sense for floating-point
// types, in which case double is fine.
constexpr double CG_SOLVER_TOL = 1.0e-14;

template<typename T, typename dataType>
Vector<dataType> cg_solve(const access_sparse_matrix<T, dataType> &A,
                          const Vector<dataType> &rhs)
{
  Vector<dataType> x0(rhs.size(), dataType(0));
  return cg_solve(A,rhs,x0);
}


template<typename T, typename dataType>
Vector<dataType> cg_solve(const access_sparse_matrix<T, dataType> &A,
                          const Vector<dataType> &rhs,
                          const Vector<dataType> &x0) {
  
  
  Vector<dataType> x(x0);
  Vector<dataType> r = rhs - A*x0;
  
  Vector<dataType> p(r);
  
  std::size_t k=0;
  std::size_t maxiter = x.size();
  
  dataType rTr_old = r.dot(r);
  double rTr_0 = rTr_old;
  
  
  if(rTr_old == dataType(0)) {
    return x; // <-- initial guess solved system exactly
  }
  
  while(k++ < maxiter) {

    auto Ap = A*p;
    dataType alpha = rTr_old/p.dot(Ap);
    x += p*alpha;
    
    r += -alpha*Ap;
    
    dataType rTr_new = r.dot(r);

    //convergence check
    if(rTr_new/rTr_0 < CG_SOLVER_TOL) {
      break;
    }

    dataType beta = rTr_new/rTr_old;
    
    rTr_old = rTr_new;
    
    p *= beta;
    p += r;
    
  } // end while
  
  return x;
}






YAFEL_NAMESPACE_CLOSE


#endif
