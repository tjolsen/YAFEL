#include "lin_alg/cg_solver.hpp"
#include "omp.h"
#include <cmath>
#include <iostream>

YAFEL_NAMESPACE_OPEN

Vector pcg_solve(const sparse_csr & A, const Vector &rhs, const Vector &x0, const Preconditioner &P) {
  
  Vector x(x0);
  Vector r = rhs - A*x0;
  Vector z = P.MinvV(r);
  Vector p(z);

  double rTr_old = r.dot(r);
  double rTz_old = r.dot(z);
  double rTr_0 = rTr_old;

  if(rTr_old == 0.0) {
    return x;
  }

  unsigned k = 0;
  unsigned maxiter = x.getLength()*2;
  while(k++ < maxiter) {
    //std::cerr << k << "/" << maxiter << std::endl;
    Vector Ap = A*p;
    double alpha = rTz_old/(p.dot(Ap));
    
    x += p*alpha;
    r += Ap*(-alpha);
    
    double rTr_new = r.dot(r);
    std::cout << k << "," << rTr_new/rTr_0 << std::endl;
    if(rTr_new/rTr_0 < CG_SOLVER_TOL) {
      break;
    }
    
    z = P.MinvV(r);

    double rTz_new = r.dot(z);
    double beta = rTz_new/rTz_old;
    p*=beta;
    p+=z;

    //reset values for next iter
    rTz_old = rTz_new;
    rTr_old = rTr_new;
  }//end while

  return x;
}


Vector cg_solve(const sparse_csr & A, const Vector & rhs) {
  Vector x(rhs.getLength(),0);
  return cg_solve(A, rhs, x);
}

Vector cg_solve(const sparse_csr & A, const Vector & rhs, const Vector &x0) {
  
  Vector x(x0);
  Vector r = rhs - A*x0;
  Vector p(r);
  
  double rTr_old = r.dot(r);
  double rTr_0 = rTr_old;
  
  if(rTr_old == 0.0) {
    return x;
  }
  
  unsigned k = 0;
  unsigned maxiter = x.getLength()*2;
  while( k++ < maxiter) {
    //std::cerr << k << "/" << maxiter << std::endl;  
    Vector Ap = A*p;
    double alpha = rTr_old/p.dot(Ap);
    x += p*alpha;
    
    r += Ap*(-alpha);
    
    double rTr_new = r.dot(r);
    
    // convergence tracking
    std::cout << k << "," << rTr_new/rTr_0 << std::endl;
    if(rTr_new/rTr_0 < CG_SOLVER_TOL) {
      //std::cerr << "converged in " << k << "steps." << std::endl;
      break;
    }
    
    double beta = rTr_new/rTr_old;
    
    rTr_old = rTr_new;
    
    p *= beta;
    p += r;
    
  }//end while
  
  return x;
}


Vector cg_solve(const FullMatrix & A, const Vector & rhs) {
  
  Vector r(rhs);
  Vector p(rhs);
  Vector x(rhs.getLength(), 0);
  
  double rTr_old = r.dot(r);
  
  unsigned k = 0;
  while( k < x.getLength()*2) {
    Vector Ap = A*p;
    double alpha = rTr_old/p.dot(Ap);
    x += p*alpha;
    
    Ap *= -alpha;
    r += Ap;
    
    double rTr_new = r.dot(r);
    
    if(rTr_new < CG_SOLVER_TOL) {
      break;
    }
    
    double beta = rTr_new/rTr_old;
    
    rTr_old = rTr_new;
    
    p *= beta;
    p += r;
    
    k++;
  }
  
  return x;
}

YAFEL_NAMESPACE_CLOSE
