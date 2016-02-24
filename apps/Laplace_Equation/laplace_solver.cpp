#include "yafel_globals.hpp"

//include mesh
#include "mesh/GenericMesh.hpp"
#include "mesh/GmshMesh.hpp"
#include "mesh/RectilinearMesh.hpp"

//include linear algebra
#include "lin_alg/Vector.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/sparse_coo.hpp"
#include "lin_alg/operators.hpp"
#include "lin_alg/solver/solvers.hpp"
#include "lin_alg/tensor/Tensor.hpp"

//include elements
#include "element/ElementFactory.hpp"
#include "element/LinQuad.hpp"
#include "element/LinTri.hpp"

//include utilities
#include "utils/DoFManager.hpp"
#include "utils/ElementType.hpp"
#include "utils/DirBC.hpp"


#include <iostream>
#include <vector>


using namespace yafel;

constexpr unsigned NSD = 2;


int main() {

  /*
   * Problem Parameters
   */
  double L = 1; // box length
  double dim_elem = 10; // elements along each dimension

  //set up basic structures
  RectilinearMesh<NSD> M(std::vector<double>(NSD,L), std::vector<std::size_t>(NSD,dim_elem));
  DoFManager dofm;
  ElementFactory<RectilinearMesh<NSD>, NSD> EF(M,dofm);

  sparse_coo<double> COO;
  Vector<double> Fsys(EF.n_dof());
  

  // assemble system
  for(std::size_t elnum=0; elnum<M.n_elements(); ++elnum) {
    std::cout << elnum << std::endl;
    Element<NSD> & el = EF.getElement(elnum);
    std::cout << "Got element\n";
    if(el.element_type == ElementType::NULL_ELEMENT) {
      continue;
    }
    std::cout << "updating...";
    el.update_element(M,elnum);
    std::cout << "done\n";
    std::size_t dof_per_el = el.dof_per_el;
    
    Matrix<double> Kloc(dof_per_el, dof_per_el, 0);
    
    for(std::size_t qpi=0; qpi<el.n_quadPoints; ++qpi) {
      std::cout << "\tqpi="<<qpi<<std::endl;
      for(std::size_t A=0; A<dof_per_el; ++A) {
        for(std::size_t B=A; B<dof_per_el; ++B) {
          Kloc(A,B) = 0;
          Kloc(B,A) = 0;
          for(std::size_t i=0; i<NSD; ++i) {
            Kloc(A,B) += el.shape_grads[qpi](A,i)*el.shape_grads[qpi](B,i)*el.JxW(qpi);
          }
          Kloc(B,A) = Kloc(A,B);
        }
      }
    }

    for(std::size_t A=0; A<dof_per_el; ++A) {
      std::size_t GA = el.global_dofs[A];
      for(std::size_t B=0; B<dof_per_el; ++B) {
        std::size_t GB = el.global_dofs[B];
        
        COO.add(GA, GB, Kloc(A,B));
      }
    }
    
  }//end element loop

  // solve system
  


  // output system

  return 0;
}
