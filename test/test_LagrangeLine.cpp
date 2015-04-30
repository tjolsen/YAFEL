#include "yafel.hpp"
#include <iostream>

using namespace yafel;

int main(int argc, char **argv) {
  
  if(argc < 2) {
    std::cout << "Supply mesh name" << std::endl;
    return 1;
  }

  Mesh M(MeshReader::gmsh_read(argv[1]));
  DoFManager dofm(1);
  ElementFactory EF(M, dofm);;
  
  Element *e = EF.getElement(0);
  e->update_element(M,0);
  
  std::cout << "Vals" << std::endl;
  for(unsigned qpi=0; qpi<e->n_quadPoints; ++qpi) {
    std::cout << "xi = " << e->quad_points[qpi](0) << ","
	      <<" w = " << e->gauss_weights[qpi] << ": ";
    
    for(unsigned node=0; node<e->nodes_per_el; ++node) {
      std::cout << e->vals[qpi](node) << " ";
    }
    std::cout << std::endl;
  }
  
  return 0;
}
