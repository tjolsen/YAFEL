#include "utils/GaussLegendreQuadrature.hpp"
#include <cmath>


YAFEL_NAMESPACE_OPEN

GaussLegendreQuadrature::GaussLegendreQuadrature(unsigned polyOrder) :
  QuadratureRule(polyOrder)
{
  
  double PI = atan2(1,1)*4;
  double TOL = 1.0e-14;

  for(unsigned xi=0; xi<polyOrder; ++xi) {
    
    double x = cos(PI*(1+xi-0.25)/(polyOrder+0.5));
    double Pold, Pn, Pnew, dP;
    
    do {
      //compute legendre polynomial
      Pold = 1; Pn = x;
      for(unsigned n=2; n<=polyOrder; ++n) {
	Pnew = ( (2*n - 1)*x*Pn - (n-1)*Pold )/n;
	
	Pold = Pn;
	Pn = Pnew;
      }
      
      //legendre polynomial derivative
      dP = (polyOrder/(x*x - 1))*(x*Pn - Pold);
      
      x -= Pn/dP;
    } while(std::abs(Pn) > TOL);
    
    nodes[xi] = x;
    //compute weight
    double w = 2.0/( (1 - x*x)*(dP*dP) );
    weights[xi] = w;
  }
  
  
}


YAFEL_NAMESPACE_CLOSE
