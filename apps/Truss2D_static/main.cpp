#include <cstdlib>
#include <iostream>
#include "Truss2D.hpp"

int main(int argc, char **argv) {
  
  if(argc < 3) {
    std::cout << "Usage: ./Truss2D <meshfilename> <outputfilename>\n";
    exit(1);
  }
  
  Truss2D T2D(argv[1], argv[2]);

  T2D.run();

  return 0;
}
