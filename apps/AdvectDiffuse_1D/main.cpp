#include "yafel.hpp"
#include <iostream>
#include <vector>
#include <string>

using namespace yafel;

int main(int argc, char **argv) {

  //parse input args --> read mesh
  if (argc < 2) {
    std::cerr << "Supply meesh name" << std::endl;
    return 1;
  }
  Mesh M(MeshReader::gmsh_read(argv[1]));
  M.reorder_rcm();
  DoFManager dofm;
  ElementFactory EF(M, dofm);
  //=====================================================

  // problem parameters
  double Cspeed = 0.1;
  double Ddiff = .01;
  double qvol = 0;
  SpatialFunction<double> volForcing([](const Vector &x){return 1.0;});

  //======================================================

  // Set up Boundary conditions
  unsigned TAG1 = 1;
  unsigned TAG2 = 2;
  unsigned COMP = 0;
  SpatialFunction<double> bcFunc1([](const Vector &x){return 0.0;});
  SpatialFunction<double> bcFunc2([](const Vector &x){return 1.0;});
  DirBC BC1(M, dofm, TAG1, COMP, bcFunc1);
  DirBC BC2(M, dofm, TAG2, COMP, bcFunc2);
  
  //======================================================
  
  // init some global data structures
  sparse_coo Kcoo;
  Vector Fsys(M.get_n_nodes(), 0.0);
  
  
  //======================================================
  
  // assemble system
  for(unsigned elnum=0; elnum<M.get_n_elems(); ++elnum) {
    Element *e = EF.getElement(elnum);
    if(e == NULL || e->n_spaceDim != 1) {
      continue;
    }
    
    e->update_element(M,elnum);

    //allocate local data structures
    unsigned Ndof = e->dof_per_el;
    FullMatrix Kloc(Ndof, Ndof, 0.0);
    Vector Floc(Ndof, 0.0);
    
    //loop over quadrature points
    for(unsigned qpi=0; qpi<e->n_quadPoints; ++qpi) {
      
      for(unsigned A=0; A<Ndof; ++A) {
	for(unsigned B=0; B<Ndof; ++B) { // <-- not symmetric, can't shortcut this loop

	  Kloc(A,B) += Cspeed*(e->vals[qpi](A))*(e->grads[qpi](B,0)) + 
	    Ddiff*(e->grads[qpi](A,0))*(e->grads[qpi](B,0))*e->JxW(qpi);

	} //end B loop
	
	Floc(A) += e->vals[qpi](A)*qvol*e->JxW(qpi);
      }//end A loop
      
    } // end qpi loop
    
    //assemble local structures into global
    for(unsigned A=0; A<Ndof; ++A) {
      unsigned AGlobal = e->global_dofs[A];
      for(unsigned B=0; B<Ndof; ++B) {
	unsigned BGlobal = e->global_dofs[B];
	
	Kcoo.add(AGlobal, BGlobal, Kloc(A,B));
      }
      
      Fsys(AGlobal) += Floc(A);
    }
    
  } //end elnum loop

  // compress coo sparse system matrix into csr
  sparse_csr Kcsr(Kcoo);
  
  //========================================================

  // apply boundary conditions
  BC1.apply(Kcsr, Fsys); //note: BC.apply() modifies the sparse_csr and Vector arguments
  BC2.apply(Kcsr, Fsys); 

  //========================================================

  // solve system (non-symmetric, use bicgstab)
  Vector u = bicgstab(Kcsr, Fsys);
  
  //=======================================================

  // write output
  VTKOutput VO;
  VTKMesh vtkm(EF);
  VTKScalarData vtku(u, VTKObject::VTKPOINTDATA, std::string("u"));
  VO.addVTKObject(&vtkm);
  VO.addVTKObject(&vtku);
  VO.write("output.vtu");
  
  return 0;
}
