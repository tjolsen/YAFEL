#include "element/LinTet.hpp"

YAFEL_NAMESPACE_OPEN

LinTet::LinTet(int dofpn) : Element(3, 1, dofpn, dofpn*4, 10, 4) {
  Vector v(n_spaceDim, 0.0);
  
  xi_0.clear();
  v(0) = 0; v(1) = 0; v(2) = 0; xi_0.push_back(v);
  v(0) = 1; v(1) = 0; v(2) = 0; xi_0.push_back(v);
  v(0) = 0; v(1) = 1; v(2) = 0; xi_0.push_back(v);
  v(0) = 0; v(1) = 0; v(2) = 1; xi_0.push_back(v);

  quad_points.clear();
  gauss_weights.clear();
  double a = 1.0/4.0;
  double b = 1.0/6.0;
  v(0) = a; v(1) = a; v(2) = a; quad_points.push_back(v); gauss_weights.push_back(b);
  
}


double LinTet::shape_value_xi(int node, const Vector &xi) const {
  
  switch(node) {
  case 0:
    return 1-xi(0)-xi(1)-xi(2);
  case 1:
    return xi(0);
  case 2:
    return xi(1);
  case 3:
    return xi(2);
  }
  
  //should never get here
  perror("LinTet::shape_value_xi called with invalid node number");
  exit(1);
  return 0;
}

double LinTet::shape_grad_xi(int node, int component, const Vector &xi) const {
  if(node==0)
    return -1;
  else
    return (node-1)==component;
}



YAFEL_NAMESPACE_CLOSE
