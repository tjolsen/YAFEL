#ifndef _YAFEL_GAUSSLEGENDREQUADRATURE_HPP
#define _YAFEL_GAUSSLEGENDREQUADRATURE_HPP

/*
 * Gauss-Legendre quadrature in multiple dimensions. The 1-D
 * quadrature is computed in the constructor, and then is extended
 * to multiple dimensions by forming tensor products.
 *
 * The "Npoints" parameter refers to the number of points in a 1D rule.
 * The total number of quadrature points is Npoints^NSD.
 */


#include "yafel_globals.hpp"
#include "quadrature/QuadratureRule.hpp"
#include <vector>
#include <cmath>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class GaussLegendreQuadrature : public QuadratureRule<NSD> {
  
public:
  using size_type = typename QuadratureRule<NSD>::size_type;
  using value_type = typename QuadratureRule<NSD>::value_type;
  using coordinate_type = typename QuadratureRule<NSD>::coordinate_type;
  
  GaussLegendreQuadrature(size_type Npoints);

private:
  void compute_1d_points(std::vector<value_type> & _nodes,
                         std::vector<value_type> & _weights,
                         size_type polyOrder);
};

/* 
 * Implementation
 *
 * (Ctor specializations in cpp file)
 */


template<unsigned NSD>
void GaussLegendreQuadrature<NSD>::compute_1d_points(std::vector<value_type> & _nodes,
                                                     std::vector<value_type> & _weights,
                                                     size_type Npoints)
{

  size_type polyOrder = Npoints;
  double PI = atan2(1,1)*4;
  double TOL = 1.0e-14;

  // Newton-Raphson solve for xi-th root of Legendre polynomial
  for(size_type xi=0; xi<polyOrder; ++xi) {

    double x = -std::cos(PI*(1+xi-0.25)/(polyOrder+0.5));
    double Pold, Pn, Pnew, dP;

    do {
      //compute legendre polynomial
      Pold = 1; Pn = x;
      for(size_type n=2; n<=polyOrder; ++n) {
        Pnew = ( (2*n - 1)*x*Pn - (n-1)*Pold )/n;

        Pold = Pn;
        Pn = Pnew;
      }

      //legendre polynomial derivative
      dP = (polyOrder/(x*x - 1))*(x*Pn - Pold);

      x -= Pn/dP;
    } while(std::abs(Pn) > TOL);

    _nodes[xi] = x;
    //compute weight
    double w = 2.0/( (1 - x*x)*(dP*dP) );
    _weights[xi] = w;
  }


}



YAFEL_NAMESPACE_CLOSE

#endif
