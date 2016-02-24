#ifndef _YAFEL_GAUSSLOBATTOQUADRATURE_HPP
#define _YAFEL_GAUSSLOBATTOQUADRATURE_HPP

#include "yafel_globals.hpp"
#include "utils/QuadratureRule.hpp"
#include <vector>
#include <cmath>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class GaussLobattoQuadrature : public QuadratureRule<NSD> {
  
public:
  using size_type = typename QuadratureRule<NSD>::size_type;
  using coordinate_type = typename QuadratureRule<NSD>::coordinate_type;

  GaussLobattoQuadrature(size_type polyOrder);
  
};

template<unsigned NSD>
GaussLobattoQuadrature<NSD>::GaussLobattoQuadrature(unsigned Npoints) :
  QuadratureRule(Npoints)
{
  
  double PI = atan2(1,1)*4;
  double TOL = 1.0e-15;
  //newton-raphson solve for quadrature points
  for(unsigned xi=0; xi<Npoints; ++xi) {
    double x = -cos((PI*xi)/(Npoints-1));
    double xold = x;
    double Pold, Pn, Pnew;

    do {
      xold = x;
      Pold = 1;
      Pn = x;
      
      //compute legendre polynomial value at current x
      for(unsigned n=2; n<Npoints; ++n) {
	Pnew = ( (2*n - 1)*x*Pn - (n-1)*Pold )/n;
	
	Pold = Pn;
	Pn = Pnew;
      }
      
      x = x - (x - Pold/Pn)/Npoints;

    } while(std::abs(x-xold)>TOL);

    nodes[xi] = x;
    //compute weight from formula
    double w = 2/(Npoints*(Npoints-1)*Pn*Pn);
    weights[xi] = w;
  }
  
}


YAFEL_NAMESPACE_CLOSE

#endif
