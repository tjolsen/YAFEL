#include "element/Element.hpp"
#include <stdio.h>

YAFEL_NAMESPACE_OPEN

Element::Element(unsigned nsd, unsigned nqp, unsigned dofpn, unsigned dofpe, int vtktype, unsigned nodespe) :
  n_spaceDim(nsd), n_quadPoints(nqp), dof_per_node(dofpn), 
  dof_per_el(dofpe), vtk_type(vtktype), nodes_per_el(nodespe)
{
  
}

FullMatrix Element::calcJ_xi(Vector xi) {
  
  FullMatrix retMat(n_spaceDim, n_spaceDim, 0.0);

  for(unsigned i=0; i<n_spaceDim; ++i) {
    for(unsigned j=0; j<n_spaceDim; ++j) {
      for(unsigned A=0; A<nodes_per_el; ++A) {
	retMat(i,j) += shape_grad_xi(A, j, xi) * nodal_coords[A](i);
      }
    }
  }
  
  return retMat;
}

void Element::calcJacobians() {
  jacobians.clear();
  detJ.clear();

  for(unsigned i=0; i<n_quadPoints; ++i) {
    jacobians.push_back( calcJ_xi(quad_points[i]) );
    
    detJ.push_back(jacobians[i].det());
  }

}

void Element::calcGrads() {
  
  grads.clear();
  
  for(unsigned qpi=0; qpi<n_quadPoints; ++qpi) {
    Vector qp = quad_points[qpi];
    
    FullMatrix Jinv = jacobians[qpi].getInverse();
    FullMatrix NG(nodes_per_el, n_spaceDim, 0.0);

    for(unsigned A=0; A<nodes_per_el; ++A) {
      for(unsigned j=0; j<n_spaceDim; ++j) {
	for(unsigned k=0; k<n_spaceDim; ++k) {
	  NG(A,j) += shape_grad_xi(A,k,qp)*Jinv(k,j);
	}
      }
    }
    
    grads.push_back(NG);
  }
  
}

void Element::calcVals() {
  vals.clear();
  for(unsigned qpi=0; qpi<n_quadPoints; ++qpi) {
    Vector qp = quad_points[qpi];
    Vector V(nodes_per_el, 0.0);
    for(unsigned A=0; A<nodes_per_el; ++A) {
      V(A) = shape_value_xi(A, qp);
    }
    vals.push_back(V);
  }
}

void Element::update_element(const Mesh &M, unsigned elnum) {

  unsigned Nnodes = M.elements[elnum].size();

  global_dofs.clear();
  nodal_coords.clear();
  element = M.elements[elnum];
  for(unsigned i=0; i<Nnodes; ++i) {
    unsigned nodeNum = M.elements[elnum][i];
    nodal_coords.push_back(M.nodal_coords[nodeNum]);
    for(unsigned j=0; j<dof_per_node; ++j) {
      unsigned dofNum = nodeNum*dof_per_node + j;
      global_dofs.push_back(dofNum);
    }
  }

  calcJacobians();
  calcGrads();
}



YAFEL_NAMESPACE_CLOSE
