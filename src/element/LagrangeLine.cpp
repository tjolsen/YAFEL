#include "element/LagrangeLine.hpp"
#include "utils/DualNumber.hpp"
#include "utils/GaussLegendreQuadrature.hpp"
#include <cmath>

YAFEL_NAMESPACE_OPEN

LagrangeLine::LagrangeLine(const DoFManager &dofm, unsigned polynomialOrder) : 
  Element(dofm, 1, polynomialOrder+1, dofm.getDofPerNode()*(polynomialOrder+1),
	  4, polynomialOrder+1)
{
  
  Vector v(n_spaceDim, -1.0);
  // lay out reference points
  double dx = 2.0/polynomialOrder;
  xi_0.clear();
  xi_0.push_back(v);
  v(0) = 1.0; xi_0.push_back(v);
  for(unsigned i=1; i<polynomialOrder; ++i) {
    v(0) = -1.0 + i*dx; xi_0.push_back(v);
  }
  
  
  GaussLegendreQuadrature GQ(polynomialOrder+1);
  
  //assign gauss points and weights
  quad_points.clear();
  gauss_weights.clear();
  
  for(unsigned i=0; i<GQ.get_nqp(); ++i) {
    v(0) = GQ.node(i);
    quad_points.push_back(v);
    gauss_weights.push_back(GQ.weight(i));
  }
  
}

double LagrangeLine::shape_value_xi(unsigned node, const Vector &xi) const {
  
  double prod = 1;
  
  for(unsigned i=0; i<nodes_per_el; ++i) {
    if(i != node)
      prod *= (xi(0)-xi_0[i](0))/(xi_0[node](0)-xi_0[i](0));
  }
  
  return prod;
}

double LagrangeLine::shape_grad_xi(unsigned node, unsigned comp, const Vector &xi) const {
  
  DualNumber<double> x(xi(comp),1.0), prod(1.0, 0.0);
  
  for(unsigned i=0; i<nodes_per_el; ++i) {
    if(i != node)
      prod = prod*(x - xi_0[i](0))/(xi_0[node](0)-xi_0[i](0));
  }
  
  return prod.second;
}

YAFEL_NAMESPACE_CLOSE
