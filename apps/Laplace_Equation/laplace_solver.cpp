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

//include output tools
#include "output/MatrixVisualization.hpp"

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

  //manually construct boundary conditions
  std::vector<std::size_t> bcnodes;
  std::vector<double> bcvals;
  std::vector<bool> bcmask(M.n_nodes(),false);

  double x0_bcval = 1;
  double x1_bcval = 0;
  double y0_bcval = 0;
  double y1_bcval = 0;

  for(std::size_t i=0; i<M.n_nodes(); ++i) {
    Tensor<NSD,1> x = M.node(i);
    if(x(0) == 0) {
      bcnodes.push_back(i);
      bcvals.push_back(x0_bcval);
      bcmask[i] = true;
    }
    else if(x(0) == L) {
      bcnodes.push_back(i);
      bcvals.push_back(x1_bcval);
      bcmask[i] = true;
    }
    else if(x(1) == 0) {
      bcnodes.push_back(i);
      bcvals.push_back(y0_bcval);
      bcmask[i] = true;
    }
    else if(x(1) == L) {
      bcnodes.push_back(i);
      bcvals.push_back(y1_bcval);
      bcmask[i] = true;
    }
  }

  DirBC BC(bcnodes, bcvals, bcmask);

  // solve system
  sparse_csr<double> Ksys(COO);
  Vector<double> Fsys(Ksys.rows(), 0);

  BC.apply(Ksys, Fsys);
  
  Vector<double> U = cg_solve(Ksys, Fsys, BC.get_ubc());

  // output system
  VTKOutput vout;

  
  
  return 0;
}
