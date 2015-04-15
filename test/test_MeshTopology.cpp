#include "yafel.hpp"
#include <iostream>
#include <string>

using namespace yafel;

int main() {
  
  Mesh M = MeshReader::gmsh_read(std::string("minsquare.msh"));
  std::cout << "Got mesh" << std::endl;
  MeshTopology MT(M);

  for(unsigned elnum=0; elnum<M.get_n_elems(); ++elnum) {
    std::cout << elnum << ": ";
    for(unsigned i=0; i<M.elements[elnum].size(); ++i) {
      unsigned neighbornum = MT.getCellNeighbor(elnum,i);
      if(neighbornum != elnum)
	std::cout << neighbornum << " ";
    }
    std::cout << std::endl;
    
  }


  return 0;
}
