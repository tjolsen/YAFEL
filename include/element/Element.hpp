#ifndef _YAFEL_ELEMENT_HPP
#define _YAFEL_ELEMENT_HPP

#include <vector>

#include "yafel_globals.hpp"
#include "lin_alg/FullMatrix.hpp"
#include "lin_alg/Vector.hpp"
#include "mesh/Mesh.hpp"

YAFEL_NAMESPACE_OPEN

class Element {
  
public:
  int n_spaceDim;
  int n_quadPoints;
  int dof_per_node;
  int dof_per_el;
  int vtk_type;
  int el_num;
  int nodes_per_el;
  
  std::vector<double> gauss_weights;
  std::vector<Vector> quad_points;
  std::vector<Vector> xi_0;
  std::vector<FullMatrix> jacobians;
  std::vector<double> detJ;
  // vector of values of dofs at qp's
  // uses Vector object to since dof's often represent a mathematical vector.
  // std::vector is used for general container.
  std::vector<Vector> vals; //what was i thinking here???
  std::vector<FullMatrix> grads; //vector of n_qp FullMatrix objects, hold grads of dof at qp's
  std::vector<int> element;
  std::vector< Vector > nodal_coords;
  std::vector<int> global_dofs;
  
  Element(int nsd, int nqp, int dofpn, int dofpe, int vtktype, int nodespe);

  // Virtual functions, specialized in child classes
  virtual ~Element();
  virtual double shape_value_xi(int node, const Vector &xi) const = 0;
  virtual double shape_grad_xi(int node, int component, const Vector &xi) const = 0;

  //Functions in base class
  FullMatrix calcJ_xi(Vector xi);
  void calcJacobians(); // calcualte Jacobians at Gauss points and store in 
  void calcGrads(); // calculate shape function gradients (wrt spatial coords) and store in "grads[qpi](A, i)"
  void update_element(const Mesh & M, int elnum);
  inline double JxW(int qpi) const { return gauss_weights[qpi]*detJ[qpi]; }
  
  
  //utility functions, might need to use in program so make public
  inline int getComp(int dof) const { return (dof % dof_per_node); }
  inline int getBase(int dof) const { return (dof/dof_per_node); }
  
  double xval(int dof, Vector xi);
  
  
};

YAFEL_NAMESPACE_CLOSE

#endif //#ifndef _YAFEL_ELEMENT_HPP
