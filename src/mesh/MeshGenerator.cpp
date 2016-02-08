#include "mesh/MeshGenerator.hpp"
#include <cstdio>
#include <cstdlib>
#include <cassert>

YAFEL_NAMESPACE_OPEN

MeshGenerator::MeshGenerator(): 
  MeshGenerator(2) {}

MeshGenerator::MeshGenerator(size_type NSD) : 
  MeshGenerator(NSD, 1.0) {}

MeshGenerator::MeshGenerator(size_type NSD, double l) : 
  MeshGenerator(NSD, l, 10) {}

MeshGenerator::MeshGenerator(size_type NSD, double l, size_type nx) :
  MeshGenerator(NSD, std::vector<double>(NSD, l), std::vector<size_type>(NSD, nx)) {}

MeshGenerator::MeshGenerator(size_type NSD, std::vector<double> l) : 
  MeshGenerator(NSD, l, std::vector<size_type>(NSD, 10)) {}

MeshGenerator::MeshGenerator(size_type NSD, std::vector<double> l, std::vector<size_type> nx) :
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
  size_type nNodes = 1;
  for(size_type i=0; i<nsd; ++i) {
    nNodes *= Nx[i];
  }
  
  std::vector<coordinate_type> nodes;
  std::vector<std::vector<size_type> > elems;
  std::vector<size_type> eltype;
  std::vector<std::vector<size_type> > tags;
  
  if(nsd == 1) {
    
    //element length
    double le = L[0]/(Nx[0]-1);
    
    //Points
    nodes.push_back(coordinate_type());
    for(size_type i=1; i<nNodes; ++i) {
      coordinate_type v;
      v(0) = i*le;
      nodes.push_back(v);
    }
    
    //Elements
    //boundary
    std::vector<size_type> e(1,0);
    elems.push_back(e);
    e[0] = nNodes-1;
    elems.push_back(e);
    tags.push_back(std::vector<size_type>(1,1));
    tags.push_back(std::vector<size_type>(1,2));
    eltype.push_back(15);
    eltype.push_back(15);
    //bulk
    for(size_type i=0; i<nNodes-1; ++i) {
      std::vector<size_type> el{i, i+1};
      elems.push_back(el);
      tags.push_back(std::vector<size_type>(1,0));
      eltype.push_back(1);
    }
    
  }
  else if(nsd == 2) {
    //Points

    double lx = L[0]/(Nx[0]-1);
    double ly = L[1]/(Nx[1]-1);
    for(size_type n=0; n<nNodes; ++n) {
      size_type r = n/Nx[0];
      size_type c = n%Nx[0];
      
      double x = c*lx;
      double y = r*ly;
      coordinate_type v;
      v(0) = x; v(1) = y;
      nodes.push_back(v);
    }

    //Elements
    //corners
    elems.push_back(std::vector<size_type>(1,0));
    tags.push_back(std::vector<size_type>{1});
    eltype.push_back(15);
    elems.push_back(std::vector<size_type>(1,Nx[0]-1));
    tags.push_back(std::vector<size_type>{2});
    eltype.push_back(15);
    elems.push_back(std::vector<size_type>(1,(Nx[1]-1)*Nx[0]));
    tags.push_back(std::vector<size_type>{3});
    eltype.push_back(15);
    elems.push_back(std::vector<size_type>(1,nNodes-1));
    tags.push_back(std::vector<size_type>{4});
    eltype.push_back(15);
    
    // edges
    for(size_type i=0; i<Nx[0]-1; ++i) {
      elems.push_back(std::vector<size_type>{i, i+1});
      tags.push_back(std::vector<size_type>{5});
      eltype.push_back(1);
    }
    for(size_type i=0; i<Nx[0]-1; ++i) {
      size_type offset = (Nx[1]-1)*Nx[0];
      elems.push_back(std::vector<size_type>{i+offset, i+1+offset});
      tags.push_back(std::vector<size_type>{6});
      eltype.push_back(1);
    }
    for(size_type i=0; i<Nx[1]-1; ++i) {
      elems.push_back(std::vector<size_type>{i*Nx[0], (i+1)*Nx[0]});
      tags.push_back(std::vector<size_type>{7});
      eltype.push_back(1);
    }
    for(size_type i=0; i<Nx[1]-1; ++i) {
      size_type offset = Nx[0]-1;
      elems.push_back(std::vector<size_type>{i*Nx[0]+offset, (i+1)*Nx[0]+offset});
      tags.push_back(std::vector<size_type>{8});
      eltype.push_back(1);
    }
    
    //bulk
    for(size_type r=0; r<Nx[1]-1; ++r) {
      for(size_type c=0; c<Nx[0]-1; ++c) {
	size_type i = r*Nx[0] + c;
	elems.push_back(std::vector<size_type>{i, i+1, i+Nx[0]+1, i+Nx[0]});
	tags.push_back(std::vector<size_type>{0});
	eltype.push_back(3);
      }
    }


  }
  else if(nsd == 3) {
    //not supporting this one just yet
    assert(false && "MeshGenerator: NSD=3 not yet supported");
    
  }
  
  
  return Mesh(nodes, elems, eltype, tags);
}
  
YAFEL_NAMESPACE_CLOSE
