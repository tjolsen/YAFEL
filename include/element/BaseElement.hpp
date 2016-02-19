#ifndef __YAFEL_BASEELEMENT_HPP
#define __YAFEL_BASEELEMENT_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/Tensor"
#include "lin_alg/Vector.hpp"
#include "lin_alg/Matrix.hpp"
#include "mesh/GenericMesh.hpp"
#include "util/DoFManager.hpp"

#include <vector>


YAFEL_NAMESPACE_OPEN

template<typename T, unsigned NSD>
class BaseElement {

public:
  using size_type = typename Tensor<NSD,1>::size_type;
  using value_type = typename Tensor<NSD,1>::value_type;
  
  /*
   * Element Interface
   */ 

  //element information
  inline size_type n_sd() const {return static_cast<T const&>(*this).n_sd();}
  inline size_type dof_per_node() const {return static_cast<T const&>(*this).dof_per_node();}

  //return vector of nodes comprising the facenum-th face of the element
  std::vector<size_type> face(size_type facenum) const {return static_cast<T const&>(*this).face(facenum);}
  size_type n_faces() const {return static_cast<T const&>(*this).n_faces();}

  
  //quadrature information
  inline size_type n_qp() const {return static_cast<T const&>(*this).n_qp();} // # of quad. pts.
  inline Tensor<NSD,1> quadrature_point(size_type qpi) const {
    return static_cast<T const&>(*this).quadrature_point(qpi);
  }
  inline value_type quadrature_weight(size_type qpi) const {
    return static_cast<T const&>(*this).quadrature_weight(qpi);
  }
  
  // shape function value and gradient at a point
  value_type shape_value_xi(size_type dof, const Tensor<NSD,1> &xi) const {
    return static_cast<T const&>(*this).shape_value_xi(dof,xi);
  }
  value_type shape_gradient_xi(size_type dof, size_type direction, const Tensor<NSD,1> &xi) const {
    return static_cast<T const&>(*this).shape_gradient_xi(dof,direction,xi);
  }

  // return "component" and "node" corresponding to element degree of freedom. Offloaded to child class.
  inline size_type get_comp(size_type dof) const {
    return static_cast<T const&>(*this).get_comp(dof);
  }
  inline size_type get_node(size_type dof) const {
    return static_cast<T const&>(*this).get_node(dof);
  }
  

  /*
   * Public BaseElement methods
   */
  //update element values
  template<typename MTYPE,unsigned MNSD>
  update_element(const GenericMesh<MTYPE,MNSD> &M, const DoFManager &dofm, size_type elnum);
  
  //return determinant of jacobian times quadrature weight at quadpoint number qpi
  value_type JxW(size_type qpi) const {return static_cast<T const&>(*this).JxW(qpi);}
  
  /*
   * CRTP glue
   */
  operator T&() {return static_cast<T&>(*this);}
  operator T const&() const {return static_cast<T const&>(*this);}
  
protected:
  //populate via construction
  size_type _dof_per_node;

  //populate via update_element
  std::vector<Vector<double> > _shape_vals; //shape function values at quadrature points [qp](dof)
  std::vector<Matrix<double> > _shape_grads; //shape function gradients at quadrature points [qp](dof, dim);
  std::vector<Tensor<NSD,2> > _jacobians; //jacobian of mapping from reference element to physical element
  std::vector<double> _detJ; //determinant of jacobian mapping
  std::vector<size_type> _element_nodes;
  std::vector<size_type> _global_dofs;


  //calculate _shape_vals
  void calc_shape_vals();

  //calculate _shape_grads
  void calc_shape_grads();

  //calculate jacobians
  void calc_jacobians();
};




/*
 * Implementation starts here
 */
template<typename T, unsigned NSD, typename MTYPE, unsigned MNSD>
BaseElement<NSD>::update_element(const GenericMesh<MTYPE,MNSD> &M,
                                 const DoFManager &dofm,
                                 size_type elnum) const
{
  _global_dofs.clear();
  _element_nodes = M.element(elnum);

  //populate _global_dofs
  for(auto n : _element_nodes) {
    for(size_type i=0; i<_dof_per_node; ++i) {
      _global_dofs.push_back( dofm.index(n, i) );
    }
  }

  calc_jacobians();  
  calc_shape_vals();
  calc_shape_grads();
}



YAFEL_NAMESPACE_CLOSE

#endif
