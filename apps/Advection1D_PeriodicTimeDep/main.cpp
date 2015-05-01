#include "yafel.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <unistd.h>

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
  double Cspeed = 10;
  double Ddiff = 0.01;
  double qvol = 0.0;
  double Tfinal = 10;
  double dt = 0.01;
  double alpha = 0.5;
  double omega = 10;
  SpatialFunction<double> volForcing([](const Vector &x)->double {return std::abs(x(0)-0.5)<.05;});

  //=====================================================
  // set up periodic boundary conditions
  unsigned LTAG = 1;
  unsigned RTAG = 2;
  unsigned Lnode = 0;
  unsigned Rnode = 0;

  for(unsigned e=0; e<M.get_n_elems(); ++e) {
    if(M.el_tags[e][0] == LTAG)
      Lnode = M.elements[e][0];
    if(M.el_tags[e][0] == RTAG)
      Rnode = M.elements[e][0];
  }
  
  //=====================================================
  
  sparse_coo Acoo, Bcoo;
  
  Vector Usol(M.get_n_nodes(), 0.0);
  Vector X(M.get_n_nodes(), 0.0);
  //=====================================================
  
  for(unsigned elnum=0; elnum<M.get_n_elems(); ++elnum) {
    Element *e = EF.getElement(elnum);
    if(e == nullptr || e->n_spaceDim != 1) {
      continue;
    }
    
    e->update_element(M,elnum);

    unsigned Ndof = e->dof_per_el;
    FullMatrix Kloc(Ndof,Ndof,0.0), Mloc(Ndof, Ndof, 0.0);
    
    for(unsigned qpi=0; qpi<e->n_quadPoints; ++qpi) {
      
      for(unsigned A=0; A<Ndof; ++A) {
	for(unsigned B=0; B<Ndof; ++B) {
	  Kloc(A,B) += -( Cspeed*(e->vals[qpi](A))*(e->grads[qpi](B,0)) 
			  + Ddiff*(e->grads[qpi](A,0))*(e->grads[qpi](B,0))
			  - (e->vals[qpi](A))*qvol
			  )*e->JxW(qpi);
	  
	  Mloc(A,B) += (e->vals[qpi](A))*(e->vals[qpi](B))*(e->JxW(qpi));
	}
      }
    }//end qpi
    
    
    //assemble into global
    for(unsigned A=0; A<Ndof; ++A) {
      unsigned GA = e->global_dofs[A];
      GA = (GA==Rnode) ? Lnode : GA;
      for(unsigned B=0; B<Ndof; ++B) {
	unsigned GB = e->global_dofs[B];
	GB = (GB==Rnode) ? Lnode : GB;
	Acoo.add(GA,GB, Mloc(A,B) - (1-alpha)*Kloc(A,B)*dt/2);
	Bcoo.add(GA,GB, Mloc(A,B) + alpha*Kloc(A,B)*dt/2);
      }
      //set initial condition for this element
      
      X(GA) = e->nodal_coords[A](0);
      Usol(GA) = std::abs(X(GA)-0.5);//sin(4*atan2(1,1)*(e->nodal_coords[A](0)));
    }
    


  }//end elnum
  
  sparse_csr Asys(Acoo), Bsys(Bcoo);

  

  MatrixVisualization MV;
  MV.spy(Bsys); return 0;
  //MV.scatter_xy(X,Usol);
  
  unsigned Ntimesteps = (unsigned)(Tfinal/dt);
  for(unsigned ti=0; ti<Ntimesteps; ++ti) {
    
    Usol = bicgstab(Asys, Bsys*Usol);
    Usol(Lnode) = Usol(Rnode);
    MV.scatter_xy(X,Usol);
    usleep(10000);
  }

  return 0;
}
