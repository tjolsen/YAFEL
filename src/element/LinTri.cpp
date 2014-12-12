#include "element/LinTri.hpp"

YAFEL_NAMESPACE_OPEN

LinTri::LinTri(unsigned dofpn) : Element(2, 1, dofpn, dofpn*3, 5, 3)
{
  
  Vector v(n_spaceDim,0.0);
  //set up parent element xi_0 vectors
  xi_0.clear();
  v(0) = 0; v(1) = 0; xi_0.push_back(v);
  v(0) = 1; v(1) = 0; xi_0.push_back(v);
  v(0) = 0; v(1) = 2; xi_0.push_back(v);

  //assign gauss points and weights
  quad_points.clear();
  gauss_weights.clear();
  /*
    //3 integration points
    double a = 1.0/6.0;
    double b = 2.0/3.0;
    v(0) = a; v(1) = a; quad_points.push_back(v); gauss_weights.push_back(a);
    v(0) = b; v(1) = a; quad_points.push_back(v); gauss_weights.push_back(a);
    v(0) = a; v(1) = b; quad_points.push_back(v); gauss_weights.push_back(a);
  */
  
  //single integration point
  double c = 1.0/3.0;
  v(0) = c; v(1) = c; quad_points.push_back(v); gauss_weights.push_back(1.0/2.0);
  
}

double LinTri::shape_value_xi(unsigned node, const Vector &xi) const {
  
  switch(node) {
  case 0: return (1 - xi(0) - xi(1));
  case 1: return xi(0);
  case 2: return xi(1);
  }
  
  return 0.0;
}

double LinTri::shape_grad_xi(unsigned node, unsigned comp, const Vector &xi) const {
  if(node==0)
    return -1;
  else
    return 1.0*((node-1)==comp);
}



YAFEL_NAMESPACE_CLOSE
