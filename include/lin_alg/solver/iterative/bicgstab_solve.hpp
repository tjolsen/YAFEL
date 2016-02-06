#ifndef __YAFEL_BICGSTAB_SOLVE_HPP
#define __YAFEL_BICGSTAB_SOLVE_HPP

/*
 * Implementation of the biconjugate gradient stabilized iterative solver for sparse matrices.
 * This algorithm involves an additional spmv over CG, but it converges for non-symmetric
 * linear systems.
 */

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/access_sparse_matrix.hpp"
#include "lin_alg/operators.hpp"

YAFEL_NAMESPACE_OPEN

constexpr double BICGSTAB_SOLVER_TOL = 1.0e-14;

template<typename T, typename dataType>
Vector<dataType> bicgstab_solve(const access_sparse_matrix<T,dataType> &A,
                                const Vector<dataType> &rhs) {
  Vector<dataType> x0(rhs.size(), dataType(0));
  return bicgstab_solve(A,rhs,x0);
}


template<typename T, typename dataType>
Vector<dataType> bicgstab_solve(const access_sparse_matrix<T,dataType> &A,
                                const Vector<dataType> &rhs,
                                const Vector<dataType> &x0) {


  Vector<dataType> r_old(rhs - A*x0),
    rhat(r_old), x(x0), v(x0.size(), dataType(0)),
    p(x0.size(), dataType(0));
    
  dataType rho_old, rho_new, alpha, w_old, w_new;
  
  std::size_t conv_check_freq = (rhs.size()/100 > 1) ? rhs.size()/100 : 1;
  
  dataType residual_0 = rhs.dot(rhs);
  
  if(residual_0 < BICGSTAB_SOLVER_TOL) {
    return rhs;
  }
  
  rho_old = dataType(1);
  alpha = dataType(1);
  w_old = dataType(1);
  
  std::size_t k=0;
  
  while(true) {

    rho_new = rhat.dot(r_old);
    
    dataType beta = (rho_new/rho_old)*(alpha/w_old);


    p = r_old + p*beta - v*(w_old*beta);

    v = A*p;    

    alpha = rho_new/rhat.dot(v);
    
    Vector<dataType> s = r_old - v*alpha;
    Vector<dataType> t = A*s;

    dataType t_dot_t = t.dot(t);
    if(t_dot_t == 0)
      w_new = 0;
    else
      w_new = t.dot(s)/t_dot_t;

    x += p*alpha + s*w_new;
    
    if(k%conv_check_freq == 0) {
      Vector<dataType> r = A*x-rhs;
      dataType res = r.dot(r);

      if(res/residual_0 < BICGSTAB_SOLVER_TOL) {
        break;
      }
    }
    
    
    r_old = s - t*w_new;
    rho_old = rho_new;
    w_old = w_new;

    ++k;
  }
  
  return x;
}


YAFEL_NAMESPACE_CLOSE


#endif
