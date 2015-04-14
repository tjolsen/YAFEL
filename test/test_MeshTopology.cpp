#include "yafel.hpp"
#include <iostream>
#include <string>

using namespace yafel;

int main() {
  
  Mesh M = MeshReader::gmsh_read(std::string("minsquare.msh"));
  std::cout << "Got mesh" << std::endl;
  MeshTopology MT(M);

  MT.print(std::cout);
  return 0;
}
