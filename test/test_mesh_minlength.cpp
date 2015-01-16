#include "yafel.hpp"
#include <iostream>
#include <string>

using namespace yafel;

int main(int argc, char **argv) {

  if(argc < 2) {
    std::cout << "Provide mesh file name.\n";
    return 1;
  }
  
  Mesh M(MeshReader::gmsh_read(std::string(argv[1])));
  std::cout << "Mesh read successful. Computing min length...\n";

  M.compute_min_length();
  
  std::cout << "Min Length = " << M.get_minLength() << std::endl;
  
  return 0;
}
