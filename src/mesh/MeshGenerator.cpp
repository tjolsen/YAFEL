#include "mesh/MeshGenerator.hpp"
#include <cstdio>
#include <cstdlib>

YAFEL_NAMESPACE_OPEN

MeshGenerator::MeshGenerator(): 
  MeshGenerator(2) {}

MeshGenerator::MeshGenerator(unsigned NSD) : 
  MeshGenerator(NSD, 1.0) {}

MeshGenerator::MeshGenerator(unsigned NSD, double l) : 
  MeshGenerator(NSD, l, 10) {}

MeshGenerator::MeshGenerator(unsigned NSD, double l, unsigned nx) :
  MeshGenerator(NSD, std::vector<double>(NSD, l), std::vector<unsigned>(NSD, nx)) {}

MeshGenerator::MeshGenerator(unsigned NSD, std::vector<double> l) : 
  MeshGenerator(NSD, l, std::vector<unsigned>(NSD, 10)) {}

MeshGenerator::MeshGenerator(unsigned NSD, std::vector<double> l, std::vector<unsigned> nx) :
  nsd(NSD), L(l), Nx(nx)
{
  if(nsd != L.size()) {
    perror("MeshGenerator(ctor): Inconsistent NSD and std::vector L length");
    exit(1);
  }
  if(nsd != Nx.size()) {
    perror("MeshGenerator(ctor): Inconsistent NSD and std::vector Nx length");
    exit(1);
  }
  if(NSD<=0 || NSD >= 4) {
    perror("MeshGenerator(ctor): Bad NSD param");
    exit(1);
  }
  
  
  
}

Mesh MeshGenerator::getMesh() {
  
  //calculate total number of nodes and elems
  unsigned nNodes = 1;
  for(unsigned i=0; i<nsd; ++i) {
    nNodes *= Nx[i];
  }
  
  std::vector<Vector> nodes;
  std::vector<std::vector<unsigned> > elems;
  std::vector<unsigned> eltype;
  std::vector<std::vector<unsigned> > tags;
  
  if(nsd == 1) {
    
    //element length
    double le = L[0]/(Nx[0]-1);
    
    //Points
    nodes.push_back(Vector(3,0.0));
    for(unsigned i=1; i<nNodes; ++i) {
      Vector v(3,0.0);
      v(0) = i*le;
      nodes.push_back(v);
    }
    
    //Elements
    //boundary
    std::vector<unsigned> e(1,0);
    elems.push_back(e);
    e[0] = nNodes-1;
    elems.push_back(e);
    tags.push_back(std::vector<unsigned>(1,1));
    tags.push_back(std::vector<unsigned>(1,2));
    eltype.push_back(15);
    eltype.push_back(15);
    //bulk
    for(unsigned i=0; i<nNodes-1; ++i) {
      std::vector<unsigned> el{i, i+1};
      elems.push_back(el);
      tags.push_back(std::vector<unsigned>(1,0));
      eltype.push_back(1);
    }
    
  }
  else if(nsd == 2) {
    //Points

    double lx = L[0]/(Nx[0]-1);
    double ly = L[1]/(Nx[1]-1);
    for(unsigned n=0; n<nNodes; ++n) {
      unsigned r = n/Nx[0];
      unsigned c = n%Nx[0];
      
      double x = c*lx;
      double y = r*ly;
      Vector v(3, 0.0);
      v(0) = x; v(1) = y;
      nodes.push_back(v);
    }

    //Elements
    //boundary
    elems.push_back(std::vector<unsigned>(1,0));
    tags.push_back(std::vector<unsigned>(1,1));
    eltype.push_back(15);
    elems.push_back(std::vector<unsigned>(1,Nx[0]-1));
    tags.push_back(std::vector<unsigned>(1,2));
    eltype.push_back(15);
    elems.push_back(std::vector<unsigned>(1,(Nx[1]-1)*Nx[0]));
    tags.push_back(std::vector<unsigned>(1,3));
    eltype.push_back(15);
    elems.push_back(std::vector<unsigned>(1,nNodes-1));
    tags.push_back(std::vector<unsigned>(1,4));
    eltype.push_back(15);
    
    for(unsigned i=0; i<Nx[0]-1; ++i) {
      elems.push_back(std::vector<unsigned>{i, i+1});
      tags.push_back(std::vector<unsigned>(1,5));
      eltype.push_back(1);
    }
    for(unsigned i=0; i<Nx[0]-1; ++i) {
      unsigned offset = (Nx[1]-1)*Nx[0];
      elems.push_back(std::vector<unsigned>{i+offset, i+1+offset});
      tags.push_back(std::vector<unsigned>(1,6));
      eltype.push_back(1);
    }
    for(unsigned i=0; i<Nx[1]-1; ++i) {
      elems.push_back(std::vector<unsigned>{i*Nx[0], (i+1)*Nx[0]});
      tags.push_back(std::vector<unsigned>(1,7));
      eltype.push_back(1);
    }
    for(unsigned i=0; i<Nx[1]-1; ++i) {
      unsigned offset = Nx[0]-1;
      elems.push_back(std::vector<unsigned>{i*Nx[0]+offset, (i+1)*Nx[0]+offset});
      tags.push_back(std::vector<unsigned>(1,8));
      eltype.push_back(1);
    }
    
    //bulk
    for(unsigned r=0; r<Nx[1]-1; ++r) {
      for(unsigned c=0; c<Nx[0]-1; ++c) {
	unsigned i = r*Nx[0] + c;
	elems.push_back(std::vector<unsigned>{i, i+1, i+Nx[0]+1, i+Nx[0]});
	tags.push_back(std::vector<unsigned>(1,0));
	eltype.push_back(3);
      }
    }


  }
  else if(nsd == 3) {
    //not supporting this one just yet
    
    
  }
  
  
  return Mesh(nodes, elems, eltype, tags);
}
  
YAFEL_NAMESPACE_CLOSE
