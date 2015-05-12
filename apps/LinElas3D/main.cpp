#include "yafel.hpp"
#include <iostream>
#include <string>
#include <cmath>

using namespace yafel;

#define LAMBDA 1.153e9
#define MU 76.9e6

inline double Cijkl(unsigned i, unsigned j, unsigned k, unsigned l) {
  return LAMBDA*(i==j && k==l) + MU*((i==k && j==l) + (i==l && j==k));
}

int main(int argc, char **argv) {
  
  if(argc<2) {
    std::cerr << "Supply mesh name\n";
    return 1;
  }

  Mesh M(MeshReader::gmsh_read(std::string(argv[1])));

  M.reorder_rcm();
  DoFManager dofm(3);
  ElementFactory EF(M,dofm);

  //=======================================================
  //Boundary conditions
  
  unsigned TAG1=1, TAG2=2;
  unsigned COMP0=0, COMP1=1, COMP2=2;
  //double alpha = 0.01;
  SpatialFunction<double> zeroFunc([](const Vector &x)->double{return 0.0;});
  SpatialFunction<double> pt1Func([](const Vector &x)->double{return 0.10;});
  SpatialFunction<double> twistFunc0([](const Vector &x)->double
				     {return -0.01*sqrt(x.dot(x))*sin(atan2(x(1),x(0)));});
  SpatialFunction<double> twistFunc1([](const Vector &x)->double
				     {return 0.01*sqrt(x.dot(x))*cos(atan2(x(1),x(0)));});
  DirBC BC1_0(M,dofm, TAG1, COMP0, zeroFunc);
  DirBC BC1_1(M,dofm, TAG1, COMP1, zeroFunc);
  DirBC BC1_2(M,dofm, TAG1, COMP2, zeroFunc);
  //DirBC BC2_0(M,dofm, TAG2, COMP0, twistFunc0);
  //DirBC BC2_1(M,dofm, TAG2, COMP1, twistFunc1);
  DirBC BC2_2(M,dofm, TAG2, COMP2, pt1Func);
  
  //=======================================================
  
  sparse_coo Kcoo;
  Vector Fsys(dofm.getNDofs(M.get_n_nodes()),0.0);
  
  
  //=======================================================  
  
  //assemble system
  for(unsigned elnum=0; elnum<M.get_n_elems(); ++elnum) {
  
    Element *e = EF.getElement(elnum);
    if(e == nullptr || e->n_spaceDim != 3) {
      continue;
    }
    
    e->update_element(M,elnum);
    
    unsigned Ndof = e->dof_per_el;
    FullMatrix Kloc(Ndof, Ndof, 0.0);
    //Vector Floc(Ndof, 0.0);
    
    for(unsigned qpi=0; qpi<e->n_quadPoints; ++qpi) {
      for(unsigned A=0; A<Ndof; ++A) {
	unsigned Abase = e->getBase(A);
	unsigned i = e->getComp(A);
	for(unsigned B=A; B<Ndof; ++B) {
	  unsigned Bbase = e->getBase(B);
	  unsigned k = e->getComp(B);
	  
	  for(unsigned j=0; j<3; ++j) {
	    for(unsigned l=0; l<3; ++l) {
	      Kloc(A,B) += e->JxW(qpi)*
		(
		 e->grads[qpi](Abase,j)*Cijkl(i,j,k,l)*e->grads[qpi](Bbase,l)
		 );
	    }
	  }
	  
	  Kloc(B,A) = Kloc(A,B);//symmetric K locally (and globally)
	}//end B
      }//end A
    }//end qpi
    
    //assemble into global
    for(unsigned A=0; A<Ndof; ++A) {
      unsigned GA = e->global_dofs[A];
      for(unsigned B=0; B<Ndof; ++B) {
	unsigned GB = e->global_dofs[B];
	Kcoo.add(GA,GB, Kloc(A,B));
      }
    }
  }//end elnum

  //compress
  sparse_csr Kcsr(Kcoo);
  //MatrixVisualization MV;
  //MV.spy(Kcsr);
  

  //=====================================================

  // apply bc's and solve
  
  BC1_0.apply(Kcsr,Fsys);
  //BC2_0.apply(Kcsr,Fsys);
  BC1_1.apply(Kcsr,Fsys);
  //BC2_1.apply(Kcsr,Fsys);
  BC1_2.apply(Kcsr,Fsys);
  BC2_2.apply(Kcsr,Fsys);
  
  Vector ubc = BC1_0.getUbc() + BC1_1.getUbc() + BC1_2.getUbc() + BC2_2.getUbc();
  //ubc += BC2_0.getUbc() + BC2_1.getUbc() + BC2_2.getUbc();
  
  ILUPreconditioner Kilu(Kcsr);
  std::cerr << "Got Kilu" << std::endl;
  //Vector u = pcg_solve(Kcsr, Fsys, ubc, Kilu);
  Vector u = cg_solve(Kcsr, Fsys, ubc);

  //=====================================================
  // write output to vtk file
  Vector v(3,0.0);
  std::vector<Vector> uOut(M.get_n_nodes(), v);
  
  for(unsigned node=0; node<M.get_n_nodes(); ++node) {
    v(0) = u(dofm.index(node,0));
    v(1) = u(dofm.index(node,1));
    v(2) = u(dofm.index(node,2));
    uOut[node] = v;
  }
  
  VTKOutput VO;
  VTKMesh vtkm(EF);
  VTKVectorData vtku(uOut, VTKObject::VTKPOINTDATA, std::string("u"));
  VO.addVTKObject(&vtkm);
  VO.addVTKObject(&vtku);
  VO.write("output.vtu");

  return 0;
}
