#ifndef __YAFEL_ELEMENT_HPP
#define __YAFEL_ELEMENT_HPP

#include <vector>

#include "yafel_globals.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include "mesh/Mesh.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class Element {
  
public:
  using size_type = typename Tensor<NSD,1>::size_type;


  size_type n_spaceDim;
  size_type n_quadPoints;
  size_type dof_per_node;
  size_type dof_per_el;
  int vtk_type;
  size_type el_num;
  size_type nodes_per_el;
  
  std::vector<double> gauss_weights;
  std::vector<Tensor<NSD,1> > quad_points;
  std::vector<Tensor<NSD,1> > xi_0;
  std::vector<Tensor<NSD,2> > jacobians;
  std::vector<double> detJ;
  std::vector<Vector<double> > vals; //vector of n_qp Vectors. holds values of shape funcs at qp's
  std::vector<Matrix<double> > grads; //vector of n_qp FullMatrix objects, hold grads of shape funcs at qp's
  std::vector<size_type> element;
  std::vector<Tensor<NSD,1> > nodal_coords;
  std::vector<size_type> global_dofs;
  DoFManager DOFM;
  
  Element(const DoFManager &dofm, size_type nsd, size_type nqp,
	  size_type dofpe, int vtktype, size_type nodespe);

  // Virtual functions, specialized in child classes
  virtual ~Element() {}
  virtual double shape_value_xi(size_type node, const Vector &xi) const = 0;
  virtual double shape_grad_xi(size_type node, size_type component, const Vector &xi) const = 0;

  //Functions in base class
  Matrix<double> calcJ_xi(const Tensor<NSD,1> &xi);
  void calcJacobians(); // calcualte Jacobians at Gauss points and store in 
  void calcGrads(); // calculate shape function gradients (wrt spatial coords) and store in "grads[qpi](A, i)"
  void calcVals(); // calculate shape function values and store in "grads[qpi](A, i)"
  void update_element(const Mesh & M, size_type elnum);
  inline double JxW(size_type qpi) const { return gauss_weights[qpi]*detJ[qpi]; }
  
  
  //utility functions, might need to use in program so make public
  inline size_type getComp(size_type dof) const { return (dof % dof_per_node); }
  inline size_type getBase(size_type dof) const { return (dof/dof_per_node); }
  
  double xval(size_type dof, const Vector & xi);
  
  
};


//================================================================
/*
 * Element Implementation starts here
 */
//================================================================




YAFEL_NAMESPACE_CLOSE

#endif //#ifndef _YAFEL_ELEMENT_HPP
