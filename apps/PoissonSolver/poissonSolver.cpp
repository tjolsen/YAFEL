#include <iostream>
#include <string>
#include "Poisson.hpp"

//using namespace yafel;

int main(int argc, char **argv) {

  if(argc < 3) {
    std::cout << "Usage: ./poissonSolver <gmesh .msh file> <outputFile.vtu>\n";
    return 1;
  }

  std::string inputFile(argv[1]);
  std::string of(argv[2]);
  
  yafel::Poisson P(argv[1]);
  P.run(of);

  return 0;
}
