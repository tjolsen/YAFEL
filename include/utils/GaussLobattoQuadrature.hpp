#ifndef _YAFEL_GAUSSLOBATTOQUADRATURE_HPP
#define _YAFEL_GAUSSLOBATTOQUADRATURE_HPP

/*
 * GaussLobattoQuadrature:
 *
 * Implementation of the Gauss-Lobatto quadrature rules.
 * In NSD>1, points were formed by tensor products of the 1D rule.
 *
 * The "Npoints" parameter refers to the number of points in a 1D rule.
 * The total number of quadrature points is Npoints^NSD.
 */

#include "yafel_globals.hpp"
#include "utils/QuadratureRule.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include <vector>
#include <cmath>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class GaussLobattoQuadrature : public QuadratureRule<NSD> {
  
public:
  using size_type = typename QuadratureRule<NSD>::size_type;
  using value_type = typename QuadratureRule<NSD>::value_type;
  using coordinate_type = typename QuadratureRule<NSD>::coordinate_type;

  GaussLobattoQuadrature(size_type polyOrder);
  
  std::vector<value_type> nodes_1d;
  std::vector<Tensor<NSD,1,size_type> > pairs;
private:
  void compute_1d_points(std::vector<value_type> &nodes, 
                         std::vector<value_type> &weights,
                         size_type Npoints);
};

/*
 * Implementation
 *
 * (Specialized ctors in cpp file)
 */

//=========================================================
template<unsigned NSD>
void GaussLobattoQuadrature<NSD>::compute_1d_points(std::vector<value_type> &_nodes, 
                                                    std::vector<value_type> &_weights,
                                                    size_type Npoints)
{
  double PI = atan2(1,1)*4;
  double TOL = 1.0e-15;
  //newton-raphson solve for quadrature points
  for(unsigned xi=0; xi<Npoints; ++xi) {
    double x = -std::cos((PI*xi)/(Npoints-1));
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
    
    _nodes[xi] = x;
    //compute weight from formula
    double w = 2/(Npoints*(Npoints-1)*Pn*Pn);
    _weights[xi] = w;
  }
  
  
}



YAFEL_NAMESPACE_CLOSE

#endif
