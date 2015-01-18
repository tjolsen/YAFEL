#include "element/LinLine.hpp"
#include <cmath>

YAFEL_NAMESPACE_OPEN

LinLine::LinLine(const DoFManager &dofm) : Element(dofm, 1, 2, 2*dofm.getDofPerNode(), 3, 2)
{
  
  Vector v(n_spaceDim, 0.0);
  //set up parent element xi_0 vectors
  xi_0.clear();
  v(0) = -1; xi_0.push_back(v);
  v(0) = 1; xi_0.push_back(v);

  //assign gauss points and weights
  quad_points.clear();
  gauss_weights.clear();

  double a = 1.0/sqrt(3.0);
  v(0) = -a; quad_points.push_back(v); gauss_weights.push_back(1.0);
  v(0) = a; quad_points.push_back(v); gauss_weights.push_back(1.0);

}

double LinLine::shape_value_xi(unsigned node, const Vector &xi) const {

  if(node == 0)
    return (1.0/2.0)*(1.0 - xi(0));
  else
    return (1.0/2.0)*(xi(0) + 1);

}

double LinLine::shape_grad_xi(unsigned node, unsigned comp, const Vector &xi) const {

  if(node == 0)
    return -1.0/2.0;
  else
    return 1.0/2.0;

}

YAFEL_NAMESPACE_CLOSE
