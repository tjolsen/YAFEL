#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "element/LinQuad.hpp"
#include "mesh/RectilinearMesh.hpp"
#include "utils/DoFManager.hpp"
#include "utils/ElementType.hpp"
#include <iostream>

using namespace yafel;

template<unsigned NSD>
bool test_1() {
  DoFManager dofm;
  LinQuad<NSD> el(dofm);
  bool good = true;

  good = good && 
    el.n_spaceDim == NSD && 
    el.n_topoDim == 2 &&
    el.n_quadPoints == 4 && 
    el.dof_per_node == 1 && 
    el.dof_per_el == 4 && 
    el.vtk_type == 9 &&
    el.element_type == ElementType::LINEAR_QUAD;

  return good;
}


//test integration over single element
template<unsigned NSD>
bool test_2() {
  
  std::vector<double> m_dims(NSD,1);
  std::vector<std::size_t> elres(NSD,1);
  
  RectilinearMesh<NSD> M(m_dims, elres);
  DoFManager dofm;

  LinQuad<NSD> el(dofm);
  el.update_element(M, 0);
  
  double sum = 0;
  // integrate f(x,y)=x over the element, should be 1/2
  // lots of things all have to be working for this to come out correctly
  for(std::size_t qpi=0; qpi<el.n_quadPoints; ++qpi) {
    auto qp = el.quad_points[qpi];
    auto x = el.xval(qp);

    for(std::size_t A=0; A<el.nodes_per_el; ++A) {
      sum += x(0)*el.shape_vals[qpi](A)*el.JxW(qpi);
    }
    
  }
  
  return abs(sum-0.5)<1.0e-8;
}

//test integration over lots of elements
template<unsigned NSD>
bool test_3() {
  
  std::vector<double> m_dims(NSD,1);
  std::vector<std::size_t> elres(NSD,10);
  
  RectilinearMesh<NSD> M(m_dims, elres);
  DoFManager dofm;

  LinQuad<NSD> el(dofm);

  
  double sum = 0;
  // integrate f(x,y)=x over the domain, should be 1/2
  // lots of things all have to be working for this to come out correctly
  for(std::size_t elnum=0; elnum<M.n_elements(); ++elnum) {
    el.update_element(M, elnum);
    for(std::size_t qpi=0; qpi<el.n_quadPoints; ++qpi) {
      auto qp = el.quad_points[qpi];
      auto x = el.xval(qp);
      
      for(std::size_t A=0; A<el.nodes_per_el; ++A) {
        sum += x(0)*el.shape_vals[qpi](A)*el.JxW(qpi);
      }
    }
  }
  return abs(sum-0.5)<1.0e-8;
}


//test gradients
template<unsigned NSD>
bool test_4() {
 
  std::vector<double> m_dims(NSD,1);
  std::vector<std::size_t> elres(NSD,1);
  
  RectilinearMesh<NSD> M(m_dims, elres);
  DoFManager dofm;

  LinQuad<NSD> el(dofm);
  
  el.update_element(M,0);
  
  bool good = true;

  good = good &&
    el.shape_grad_xi(0, 0, Tensor<NSD,1>{-1, -1}) == -1.0/2.0 &&
    el.shape_grad_xi(0, 1, Tensor<NSD,1>{-1, -1}) == -1.0/2.0 &&
    el.shape_grad_xi(2, 0, Tensor<NSD,1>{1, 1}) == 1.0/2.0 &&
    el.shape_grad_xi(2, 1, Tensor<NSD,1>{1, 1}) == 1.0/2.0;

  return good;
}

int main() {

  int retval = 0;
  if(!test_1<2>()) {
    std::cerr << "Failed test_1<2>" << std::endl;
    retval |= 1<<0;
  }
  if(!test_1<3>()) {
    std::cerr << "Failed test_1<3>" << std::endl;
    retval |= 1<<0;
  }
  if(!test_2<2>()) {
    std::cerr << "Failed test_2<2>" << std::endl;
    retval |= 1<<1;
  }
  if(!test_3<2>()) {
    std::cerr << "Failed test_3<2>" << std::endl;
    retval |= 1<<2;
  }
  if(!test_3<3>()) {
    std::cerr << "Failed test_3<3>" << std::endl;
    retval |= 1<<2;
  }
  if(!test_4<2>()) {
    std::cerr << "Failed test_4<2>" << std::endl;
    retval |= 1<<3;
  }
  
  
  
  return retval;
}
