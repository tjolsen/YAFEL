#include <iostream>
#include "Truss2D.hpp"


Truss2D::Truss2D(const char *Mfname, const char *outFname) :
MeshFilename(Mfname), OutputFilename(outFname)
{

  //material parameters for steel, 1cm x 1cm cross-section
  Eyoungs = 2.0e11; // 200 GPa
  Axsection = 1.0e-4; //m^2
  rho = 8000; //density kg/m^3
}

void Truss2D::setup() {
  
  //read mesh
  M = MeshReader::gmsh_read(MeshFilename);

  //boundary conditions
  bcnodes.clear();
  bccomps.clear();
  bcvals.clear();
  
  for(unsigned i=0; i<M.el_tags.size(); ++i) {
    
    if(M.el_tags[i][0] == 1) {
      for(unsigned j=0; j<M.elements[i].size(); ++j) {
	bcnodes.push_back(M.elements[i][j]);
	bccomps.push_back(0);
	bcvals.push_back(0.0);

	bcnodes.push_back(M.elements[i][j]);
	bccomps.push_back(1);
	bcvals.push_back(0.0);
      }
    }
    if(M.el_tags[i][0] == 2) {
      for(unsigned j=0; j<M.elements[i].size(); ++j) {
	bcnodes.push_back(M.elements[i][j]);
	bccomps.push_back(0);
	bcvals.push_back(0);

	bcnodes.push_back(M.elements[i][j]);
	bccomps.push_back(1);
	bcvals.push_back(0.0);
      }
    }
    
  }

}

void Truss2D::assemble() {
  
  unsigned nsd = 2;
  unsigned dof = 4;
  sparse_coo Ksys_coo;
  Vector F(M.get_n_nodes()*nsd, 0.0);
  
  double Grav = 9.81;
  
  for(unsigned elnum = 0; elnum < M.elements.size(); ++elnum) {
    
    if(M.element_type[elnum] != 1) { //skip non-line elements
      continue;
    }    
    
    Vector Cdir(nsd, 0.0);
    
    unsigned n1, n2;
    n1 = M.elements[elnum][0];
    n2 = M.elements[elnum][1];
    
    Vector dx = M.nodal_coords[n2] - M.nodal_coords[n1];

    double L = sqrt(dx.dot(dx));
    for(unsigned i=0; i<nsd; ++i) {
      Cdir(i) = dx(i)/L;
    }
    
    FullMatrix Kloc(dof, dof, 0.0);
    Vector LoadVec(dof,0.0);
    for(unsigned A=0; A<dof; ++A) {
      unsigned abase = A/2;
      unsigned acomp = A%nsd;
      for(unsigned B=0; B<dof; ++B) {
	unsigned bbase = B/2;
	unsigned bcomp = B%nsd;
	int factor = (abase==bbase) ? 1 : -1;

	Kloc(A,B) += factor*Cdir(acomp)*Cdir(bcomp)*Eyoungs*Axsection/L;
      }
      if(acomp == 1) {
	LoadVec(A) += -Grav*rho*L*Axsection/2;
      }
    }

    std::vector<unsigned> globaldof;
    globaldof.clear();
    globaldof.push_back(n1*nsd);
    globaldof.push_back(n1*nsd + 1);
    globaldof.push_back(n2*nsd);
    globaldof.push_back(n2*nsd + 1);
    
    //assemble local into global
    for(unsigned A=0; A<dof; ++A) {
      for(unsigned B=0; B<dof; ++B) {
	Ksys_coo.add(globaldof[A], globaldof[B], Kloc(A,B));
      }
      F(globaldof[A]) += LoadVec(A);
    }
    
  }

  Kcoo = Ksys_coo;
  Fsys = F;

}

void Truss2D::solve() {

  sparse_csr Ksys_csr(Kcoo);
  sparse_csr Ksys_nobc(Ksys_csr);
  Vector Ubc(2*M.get_n_nodes(),0.0);

  //apply BCs
  for(unsigned i=0; i<bcvals.size(); ++i) {
    unsigned dof = 2*bcnodes[i] + bccomps[i];

    Ubc(dof) = bcvals[i];
    
    Ksys_csr.zero_row(dof); //cheap
    Ksys_csr.zero_col(dof); //expensive
    Ksys_csr.assign(dof,dof,1.0); //potentially expensive

  }

  //Ksys_csr.print_sparse();


  Vector rhs = Fsys - Ksys_nobc*Ubc;  
  for(unsigned i=0; i<bcvals.size(); ++i) {
    unsigned dof = 2*bcnodes[i] + bccomps[i];
    rhs(dof) = bcvals[i];
  }

  Vector u = cg_solve(Ksys_csr, rhs, Ubc);
  Usol = u;
}

void Truss2D::postprocess() {
  
  stress.clear();
  
  for(unsigned elnum=0; elnum < M.get_n_elems(); ++elnum) {
    if(M.element_type[elnum] != 1)  {
      continue;
    }

    unsigned n1 = M.elements[elnum][0];
    unsigned n2 = M.elements[elnum][1];
    Vector u1(3,0.0);
    u1(0) = Usol(2*n1);
    u1(1) = Usol(2*n1 + 1);
    Vector u2(3,0.0);
    u2(0) = Usol(2*n2);
    u2(1) = Usol(2*n2 + 1);
    
    Vector dx = M.nodal_coords[n1] - M.nodal_coords[n2];
    double L0 = dx.dot(dx);
    
    Vector dx_curr = dx + u1 - u2;
    double L = dx_curr.dot(dx_curr);
    
    double strain = (L - L0)/L0;
    
    stress.push_back(Eyoungs * strain / 1.0e6); //MPa
    
  }

}

void Truss2D::output() {
  
  std::vector<Vector> uout;
  uout.clear();
  
  for(unsigned i=0; i<Usol.getLength(); i+=2) {
    Vector v(2,0.0);
    v(0) = Usol(i);
    v(1) = Usol(i+1);
    uout.push_back(v);
  }
  
  ElementFactory EF(M, 2);
  
  VTKOutput VO;
  VTKVectorData vtku(uout, VTKObject::VTKPOINTDATA, std::string("U"));
  VTKMesh vtkm(EF);
  VTKScalarData vtkstress(stress, VTKObject::VTKCELLDATA, std::string("stress"));
  
  VO.addVTKObject(&vtkm);
  VO.addVTKObject(&vtku);
  VO.addVTKObject(&vtkstress);
  VO.write(OutputFilename);

}


void Truss2D::run() {

  setup();
  assemble();
  solve();
  postprocess();
  output();

}
