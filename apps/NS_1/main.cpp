#include <iostream>
#include <cstdlib>

#include "yafel.hpp"
#include "ChorinMethod.hpp"

using namespace yafel;

int main(int argc, char **argv) {
  
  if(argc < 3) {
    std::cout << "Supply mesh filename and number of spatial dimensions" << std::endl;
    exit(1);
  }
  
  double Tfinal = 10;
  
  ChorinMethod CM(argv[1], "ChorinOutput", atoi(argv[2]), Tfinal);
  CM.run();
  
  return 0;
}

