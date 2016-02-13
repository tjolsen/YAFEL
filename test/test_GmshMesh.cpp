#include "yafel_globals.hpp"
#include "mesh/GmshMesh.hpp"
#include <iostream>

using namespace yafel;

bool test_1() {
  GmshMesh<3> M("square.msh");


  return true;
}


int main() {

  int retval(0);

  if(!test_1()) {
    std::cerr << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }

  return retval;
}
