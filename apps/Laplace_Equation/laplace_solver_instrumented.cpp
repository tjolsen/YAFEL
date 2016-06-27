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
#include "utils/CG_DoFManager.hpp"
#include "utils/ElementType.hpp"
#include "utils/DirBC.hpp"

//include output tools
#include "output/MatrixVisualization.hpp"
#include "output/VTKOutput.hpp"
#include "output/VTKMesh.hpp"
#include "output/VTKScalarData.hpp"

#include <iostream>
#include <vector>
#include <chrono>

using namespace yafel;

constexpr unsigned NSD = 2;


int main() {

  /*
   * Problem Parameters
   */
  double L = 1; // box length
  double dim_elem = 200; // elements along each dimension

  //set up basic structures
  RectilinearMesh<NSD> M(std::vector<double>(NSD,L), std::vector<std::size_t>(NSD,dim_elem));
  CG_DoFManager<RectilinearMesh<NSD>,NSD> dofm(M,1);
  ElementFactory<RectilinearMesh<NSD>, NSD> EF(M,dofm);

  sparse_coo<double> COO;

  // assemble system
  std::cout << "Assembling..." << std::endl;
  for(std::size_t elnum=0; elnum<M.n_elements(); ++elnum) {
    Element<NSD> & el = EF.getElement(elnum);

    if(el.element_type == ElementType::NULL_ELEMENT) {
      continue;
    }

    //update element shape values, gradients
    el.update_element(M,elnum);
    
    std::size_t dof_per_el = el.dof_per_el;

    Matrix<double> Kloc(dof_per_el, dof_per_el, 0);
    
    for(std::size_t qpi=0; qpi<el.n_quadPoints; ++qpi) {
      for(std::size_t A=0; A<dof_per_el; ++A) {
        for(std::size_t B=0; B<dof_per_el; ++B) {
          for(std::size_t i=0; i<NSD; ++i) {
            Kloc(A,B) += el.shape_grads[qpi](A,i)*el.shape_grads[qpi](B,i)*el.JxW(qpi);
          }
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
  
  std::cout << "done" << std::endl
            << "Boundary conditions..." << std::endl;
  
  //manually construct boundary conditions.
  // this will be abstracted away into something nicer at some point
  std::vector<std::size_t> bcnodes;
  std::vector<double> bcvals;
  std::vector<bool> bcmask(M.n_nodes(),false);

  std::vector<double> x0_bcvals(NSD,0), x1_bcvals(NSD,0);
  for(std::size_t i=0; i<NSD; ++i) {
    x0_bcvals[i] = i+1;
  }

  for(std::size_t i=0; i<M.n_nodes(); ++i) {
    Tensor<NSD,1> x = M.node(i);
    for(std::size_t dim=0; dim<NSD; ++dim) {

      if(x(dim) == 0) {
        bcnodes.push_back(i);
        bcvals.push_back(x0_bcvals[dim]);
        bcmask[i] = true;
      }
      else if(x(dim) == L) {
        bcnodes.push_back(i);
        bcvals.push_back(x1_bcvals[dim]);
        bcmask[i] = true;
      }
    }
  }
  
  //make dirichlet bc object using vectors created above
  DirBC BC(bcnodes, bcvals, bcmask);

  
  //set up sparse matrix data structure for solving
  sparse_csr<double> Ksys(COO);
  Vector<double> Fsys(Ksys.rows(), 0);

  //apply bc to linear system
  BC.apply(Ksys, Fsys);

  std::cout << "done" << std::endl; // bc done

  //solve linear system using Preconditioned Conjugate Gradient, with the boundary conditions
  // as an initial guess (rather than just zero)

  std::cout << "Making preconditioner..." << std::endl;
  auto t1 = std::chrono::steady_clock::now();
  ILUPreconditioner<double> precond(Ksys);
  auto t2 = std::chrono::steady_clock::now();
  double precondtime = ((double)((t2-t1).count()))*std::chrono::steady_clock::period::num/std::chrono::steady_clock::period::den * 1000;
  std::cout << "done: Constructed in " << precondtime << " ms\n";
  
  std::cout << "Solving..." << std::endl;

  double TOL = 1.0e-10; // relative residual tolerance for iterative solver

  int Nsolves = 10;

  std::vector<double> cg_times, pcg_times;
  for(int i=0; i<Nsolves; ++i) {
      std::cout << i << " "; std::cout.flush();
      auto t1 = std::chrono::steady_clock::now();

      Vector<double> U = pcg_solve(Ksys, Fsys, precond, BC.get_ubc(), TOL);

      auto t2 = std::chrono::steady_clock::now();

      Vector<double> U2 = cg_solve(Ksys, Fsys, BC.get_ubc(), TOL);

      auto t3 = std::chrono::steady_clock::now();


      auto cgdur = t3-t2;
      auto pcgdur = t2-t1;

      double cgtime = ((double)cgdur.count())*std::chrono::steady_clock::period::num/std::chrono::steady_clock::period::den * 1000;
      double pcgtime = ((double)pcgdur.count())*std::chrono::steady_clock::period::num/std::chrono::steady_clock::period::den * 1000;

      cg_times.push_back(cgtime);
      pcg_times.push_back(pcgtime);


      if(i == 0) {
	  // output system
	  VTKOutput vout;
	  VTKMesh<RectilinearMesh<NSD>, NSD> vtkm(M);
	  VTKScalarData vtku(U, VTKObject::VTKPOINTDATA, "U");
	  VTKScalarData vtkubc(BC.get_ubc(), VTKObject::VTKPOINTDATA, "U BC");
	  
	  vout.addVTKObject(&vtkm);
	  vout.addVTKObject(&vtku);
	  vout.addVTKObject(&vtkubc);
	  
	  vout.write("output.vtu");
      }


  }
  double cgavg(0), pcgavg(0);
  for(auto t : cg_times)
      cgavg += static_cast<double>(t);
  for(auto t : pcg_times)
      pcgavg += static_cast<double>(t);

  cgavg /= cg_times.size();
  pcgavg /= pcg_times.size();

  std::cout << "done" << std::endl;

  std::cout << "CG average time: " << cgavg << std::endl;
  std::cout << "PCG average time: " << pcgavg << std::endl;

  return 0;
}
