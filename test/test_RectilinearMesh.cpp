#include "yafel_globals.hpp"
#include "mesh/RectilinearMesh.hpp"
#include <iostream>

using namespace yafel;

// test n_nodes, n_elements for 1D RectlinearMesh
bool test_1() {
  RectilinearMesh<1> M;
  return (M.n_nodes()==11) && (M.n_elements()==10);
}

// test n_nodes, n_elements for 2D RectlinearMesh
bool test_2() {
  RectilinearMesh<2> M;
  return (M.n_nodes()==121) && (M.n_elements()==100);
}

// test n_nodes, n_elements for 3D RectlinearMesh
bool test_3() {
  RectilinearMesh<3> M;
  return (M.n_nodes()==1331) && (M.n_elements()==1000);
}

// test node(n) for 2D RectilinearMesh
bool test_4() {
  
  //default mesh with 10 elements in each direction, 1x1 square
  RectilinearMesh<2> M;
  bool good = true;
  Tensor<2,1,double> x;
  
  // node at (0,0)
  auto mx = M.node(0);
  good = good && contract<1>(x-mx,x-mx)==0;

  // node at (1,0)
  mx = M.node(10);
  x(0) = 1;
  good = good && contract<1>(x-mx, x-mx)==0;

  // node at (.5,5)
  mx = M.node(55);
  x(0) = .5; x(1) = .5;
  good = good && contract<1>(x-mx, x-mx)==0;

  return good;
}



int main() {
  int retval = 0;

  if(!test_1()) {
    std::cerr << "Failed test_1" << std::endl;
    retval |= 1<<0;
  }
  if(!test_2()) {
    std::cerr << "Failed test_2" << std::endl;
    retval |= 1<<1;
  }
  if(!test_3()) {
    std::cerr << "Failed test_3" << std::endl;
    retval |= 1<<2;
  }


  return retval;
}
