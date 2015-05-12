#include "lin_alg/bicgstab.hpp"
#include <cmath>
#include <iostream>

YAFEL_NAMESPACE_OPEN

Vector bicgstab(const sparse_csr & A, const Vector & rhs) {
  
  Vector r_old(rhs), 
    rhat(rhs), 
    x(rhs.getLength(), 0.0),
    v(rhs.getLength(),0.0), 
    p(rhs.getLength(), 0.0);
  
  double rho_old, rho_new, alpha, w_old, w_new;
  
  unsigned conv_check_freq = std::max((unsigned)1, rhs.getLength()/100);
  
  double residual_0 = rhs.dot(rhs);

  if(residual_0 < BICGSTAB_TOL)
    return rhs;


  rho_old = 1.0;
  alpha= 1.0;
  w_old = 1.0;
  
  unsigned k=0;
  while(true) {
    rho_new = rhat.dot(r_old);
    
    double beta = (rho_new/rho_old)*(alpha/w_old);
    
    p = r_old + p*beta - v*(w_old*beta);
    
    v = A*p;
    
    
    alpha = rho_new/rhat.dot(v);
    
    Vector s = r_old - v*alpha;
    Vector t = A*s;

    double t_dot_t = t.dot(t);
    if(t_dot_t == 0)
      w_new = 0;
    else
      w_new = t.dot(s) / t_dot_t;

    x += p*alpha + s*w_new;

    if(k % conv_check_freq == 0) {
      Vector r = A*x - rhs;
      double res = r.dot(r);

      //convergence reporting
      //std::cout << k << "," << res/residual_0 << std::endl;
      if(res/residual_0 < BICGSTAB_TOL) {
	//std::cerr << "converged in " << k << "steps.\n";
	break;
      }
    }
    
    //reset old/new vars for next iteration
    r_old = s - t*w_new;
    rho_old = rho_new;
    w_old = w_new;
    
    ++k;
  }
  
  return x;
}

YAFEL_NAMESPACE_CLOSE
