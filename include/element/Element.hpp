#ifndef _YAFEL_ELEMENT_HPP
#define _YAFEL_ELEMENT_HPP

#include <vector>

#include "yafel_globals.hpp"
#include "lin_alg/FullMatrix.hpp"
#include "lin_alg/Vector.hpp"
#include "mesh/Mesh.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

class Element {
  
public:
  unsigned n_spaceDim;
  unsigned n_quadPoints;
  unsigned dof_per_node;
  unsigned dof_per_el;
  int vtk_type;
  unsigned el_num;
  unsigned nodes_per_el;
  
  std::vector<double> gauss_weights;
  std::vector<Vector> quad_points;
  std::vector<Vector> xi_0;
  std::vector<FullMatrix> jacobians;
  std::vector<double> detJ;
  // vector of values of dofs at qp's
  // uses Vector object to since dof's often represent a mathematical vector.
  // std::vector is used for general container.
  std::vector<Vector> vals; //vector of n_qp Vectors. holds values of shape funcs at qp's
  std::vector<FullMatrix> grads; //vector of n_qp FullMatrix objects, hold grads of shape funcs at qp's
  std::vector<unsigned> element;
  std::vector< Vector > nodal_coords;
  std::vector<unsigned> global_dofs;
  DoFManager DOFM;
  
  Element(const DoFManager &dofm, unsigned nsd, unsigned nqp,
	  unsigned dofpe, int vtktype, unsigned nodespe);

  // Virtual functions, specialized in child classes
  virtual ~Element() {}
  virtual double shape_value_xi(unsigned node, const Vector &xi) const = 0;
  virtual double shape_grad_xi(unsigned node, unsigned component, const Vector &xi) const = 0;

  //Functions in base class
  FullMatrix calcJ_xi(const Vector &xi);
  void calcJacobians(); // calcualte Jacobians at Gauss points and store in 
  void calcGrads(); // calculate shape function gradients (wrt spatial coords) and store in "grads[qpi](A, i)"
  void calcVals(); // calculate shape function values and store in "grads[qpi](A, i)"
  void update_element(const Mesh & M, unsigned elnum);
  inline double JxW(unsigned qpi) const { return gauss_weights[qpi]*detJ[qpi]; }
  
  
  //utility functions, might need to use in program so make public
  inline unsigned getComp(unsigned dof) const { return (dof % dof_per_node); }
  inline unsigned getBase(unsigned dof) const { return (dof/dof_per_node); }
  
  double xval(unsigned dof, Vector xi);
  
  
};

YAFEL_NAMESPACE_CLOSE

#endif //#ifndef _YAFEL_ELEMENT_HPP
