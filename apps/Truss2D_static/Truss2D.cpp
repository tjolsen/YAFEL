#include "Truss2D.hpp"


Truss2D::Truss2D(const char *Mfname, const char *outFname) :
MeshFilename(Mfname), OutputFilename(outFname), 
{
  Eyoungs = 200e9; // 200 GPa
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
  
  for(int i=0; i<M.el_tags.size(); ++i) {
    
    if(M.el_tags[i][0] == 1) {
      for(int j=0; j<M.elements[i].size(); ++j) {
	bcnodes.push_back(M.elements[i][j]);
	bccomps.push_back(0);
	bcvals.push_back(0.0);

	bcnodes.push_back(M.elements[i][j]);
	bccomps.push_back(1);
	bcvals.push_back(0.0);
      }
    }
    if(M.el_tags[i][0] == 2) {
      for(int j=0; j<M.elements[i].size(); ++j) {
	bcnodes.push_back(M.elements[i][j]);
	bccomps.push_back(1);
	bcvals.push_back(0.0);
      }
    }
    
  }
  
}

void Truss2D::assemble() {
  
  sparse_coo Ksys_coo;
  unsigned nsd = 2;
  unsigned dof = 4;
  
  for(unsigned elnum = 0; elnum < M.elements.size(); ++elnum) {
    
    if(M.el_tags[elnum][0] != 2) //skip non-line elements
      continue;
    
    
    Vector Cdir(nsd, 0.0);
    
    unsigned n1, n2;
    n1 = M.elements[elnum][0];
    n2 = M.elements[elnum][1];
    
    Vector dx = M.nodal_coords[n2] - M.nodal_coords[n1];
    
    L = sqrt(dx.dot(dx));
    for(int i=0; i<nsd; ++i) {
      Cdir(i) = dx(i)/L;
    }

    FullMatrix Kloc(dof, dof, 0.0);
    for(unsigned A=0; A<dof; ++A) {
      unsigned abase = A/2;
      for(unsigned B=0; B<dof; ++B) {
	unsigned bbase = B/2;

	int factor = (abase==bbase) ? 1 : -1;

	Kloc(A,B) += 
      }
    }
  }
  
}
