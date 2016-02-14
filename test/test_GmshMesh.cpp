#include "yafel_globals.hpp"
#include "mesh/GmshMesh.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include <iostream>

using namespace yafel;


// read in a mesh called "square.msh" (included in repo). This is a
// 1x1 square with four triangular elements.
bool test_1() {
  GmshMesh<3> M("square.msh");
  bool good = true;

  Tensor<3,1> x0, x1, x2, x3, x4;
  x1(0) = 1;
  x2(0) = 1; x2(1) = 1;
  x3(1) = 1;
  x4(0) = 0.5; x4(1) = 0.5;

  // test correct node read
  Tensor<3,1> dx;
  dx = M.node(0) - x0;
  good = good && contract<1>(dx,dx) < 1e-10;
  dx = M.node(1) - x1;
  good = good && contract<1>(dx,dx) < 1e-10;
  dx = M.node(2) - x2;
  good = good && contract<1>(dx,dx) < 1e-10;
  dx = M.node(3) - x3;
  good = good && contract<1>(dx,dx) < 1e-10;
  dx = M.node(4) - x4;
  good = good && contract<1>(dx,dx) < 1e-10;

  std::vector<std::size_t> e0{3,0,4}, e1{1,4,0}, e2{1,2,4}, e3{2,3,4};
  
  std::vector<std::vector<std::size_t>> elems = {e0,e1,e2,e3};
  for(std::size_t i=0; i<elems.size(); ++i) {
    good = good && M._element_type[i] == 2;
    for(std::size_t j=0; j<elems[i].size(); ++j) {
      good = good && elems[i][j] == M.element(i)[j];
    }
  }
  
  return good;
}


int main() {

  int retval(0);

  if(!test_1()) {
    std::cerr << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }

  return retval;
}
