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


using namespace yafel;

constexpr unsigned NSD = 3;


int main() {

    /*
     * Problem Parameters
     */
    //double L = 1; // box length
    //double dim_elem = 120; // elements along each dimension

    //set up basic structures
    //RectilinearMesh<NSD> M(std::vector<double>(NSD,L), std::vector<std::size_t>(NSD,dim_elem));
    //CG_DoFManager<RectilinearMesh<NSD>,NSD> dofm(M,1);
    //ElementFactory<RectilinearMesh<NSD>, NSD> EF(M,dofm);
    
    GmshMesh<NSD> M("square.msh");

    CG_DoFManager<GmshMesh<NSD>,NSD> dofm(M,1);
    ElementFactory<GmshMesh<NSD>, NSD> EF(M,dofm);

    sparse_coo<double> COO;

    // assemble system
    std::cout << "Assembling..." << std::endl;
    for(std::size_t elnum=0; elnum<M.n_elements(); ++elnum) {
        Element<NSD> & el = EF.getElement(elnum);

        if(el.element_type == ElementType::NULL_ELEMENT || el.n_topoDim != NSD) {
	    continue;
        }
	
        //update element shape values, gradients
        el.update_element(M,elnum);

        std::size_t dof_per_el = el.dof_per_el;
	
        Matrix<double> Kloc(dof_per_el, dof_per_el, 0);

        //loop over quadrature points
        for(std::size_t qpi=0; qpi<el.n_quadPoints; ++qpi) {
            double jxw = el.JxW(qpi);
        
            // add to local K_e (Kloc) matrix
            for(std::size_t A=0; A<dof_per_el; ++A) {
                for(std::size_t B=0; B<dof_per_el; ++B) {

                    //dot product of shape func gradients in integrand
                    for(std::size_t i=0; i<NSD; ++i) {
                        Kloc(A,B) += el.shape_grads[qpi](A,i)*el.shape_grads[qpi](B,i)*jxw;
                    }

                }
            }
        }

        //assemble Kloc into global K matrix
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



    //Set up BCs using physical ID of boundary
    int pid_1 = 10;
    int pid_2 = 11;

    SpatialFunction<NSD,double> zero_func{[](const Tensor<NSD,1,double>&){return 0.0;} };
    SpatialFunction<NSD,double> one_func{[](const Tensor<NSD,1,double>&){return 1.0;} };
    DirBC BC(M, dofm, pid_1, 0, zero_func);
    DirBC BC2(M, dofm, pid_2, 0, one_func);
  
    //set up sparse matrix data structure for solving
    sparse_csr<double> Ksys(COO);
    Vector<double> Fsys(Ksys.rows(), 0);

    //apply bc to linear system
    BC.apply(Ksys, Fsys);
    BC2.apply(Ksys, Fsys);

    auto trips = Ksys.copy_triplets();
    for(auto t : trips) {
	std::cout << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t) << std::endl;
    }

    std::cout << "done" << std::endl; // bc done

    //solve linear system using the Conjugate Gradient method, with the boundary conditions
    // as an initial guess (rather than just zero)
    std::cout << "Solving..." << std::endl;

    double TOL = 1.0e-14; // relative residual tolerance for iterative solver
    Vector<double> U = cg_solve(Ksys, Fsys, BC.get_ubc(), TOL);

    std::cout << "done" << std::endl;

    // output system
    VTKOutput vout;
    //VTKMesh<RectilinearMesh<NSD>, NSD> vtkm(M);
    VTKMesh<GmshMesh<NSD>, NSD> vtkm(M);
    VTKScalarData vtku(U, VTKObject::VTKPOINTDATA, "U");
    VTKScalarData vtkubc(BC.get_ubc(), VTKObject::VTKPOINTDATA, "U BC");

    vout.addVTKObject(&vtkm);
    vout.addVTKObject(&vtku);
    vout.addVTKObject(&vtkubc);
  
    vout.write("output.vtu");
    return 0;
}
