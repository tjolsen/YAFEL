#include "yafel.hpp"
#include <iostream>
#include <string>

using namespace yafel;


int main(int argc, char **argv) {

  Mesh M;
  
  if(argc < 2) {
    MeshGenerator MG(2, 1.0, 10);
    M = MG.getMesh();
  }
  else {
    M = MeshReader::gmsh_read(std::string(argv[1]));
  }
  
  M.reorder_rcm();
  
  sparse_coo Acoo;
  for(auto e = M.elements.begin(); e<M.elements.end(); ++e) {
    
    for(auto i=e->begin(); i<e->end(); ++i) {
      for(auto j=e->begin(); j<e->end(); ++j) {
	Acoo.add(*i, *j, 1);
      }
    }
  }

  sparse_csr Acsr(Acoo);

  MatrixVisualization::spy(Acoo);
  
  return 0;
}
