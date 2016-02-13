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
  good = good && contract<1>(x-mx,x-mx)<1e-10;
  
  // node at (1,0)
  mx = M.node(10);
  x(0) = 1;
  good = good && contract<1>(x-mx, x-mx)<1e-10;

  // node at (.5,5)
  mx = M.node(60);
  x(0) = .5; x(1) = .5;
  good = good && contract<1>(x-mx, x-mx)<1e-10;

  return good;
}


//test element(e) for 1D RectilinearMesh
bool test_5() {
  RectilinearMesh<1> M;
  bool good = true;
  
  for(std::size_t e=0; e<M.n_elements(); ++e) {
    auto el = M.element(e);
    good = good && (el[0]==e && el[1]==e+1);
  }

  return good;
}


//test element(e) for 2D RectilinearMesh
bool test_6() {
  RectilinearMesh<2> M;
  bool good = true;

  // get first element
  auto el = M.element(0);
  std::vector<std::size_t> correct{0, 1, 12, 11};
  for(std::size_t i=0; i<4; ++i) {
    good = good && el[i]==correct[i];
  }
  
  // second element
  el = M.element(1);
  correct = std::vector<std::size_t>{1,2,13,12};
  for(std::size_t i=0; i<4; ++i) {
    good = good && el[i]==correct[i];
  }

  //element from second row
  el = M.element(15);
  correct = std::vector<std::size_t>{16,17,28,27};
  for(std::size_t i=0; i<4; ++i) {
    good = good && el[i]==correct[i];
  }


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
 if(!test_4()) {
    std::cerr << "Failed test_4" << std::endl;
    retval |= 1<<3;
  }
 if(!test_5()) {
    std::cerr << "Failed test_5" << std::endl;
    retval |= 1<<4;
  }
 if(!test_6()) {
    std::cerr << "Failed test_6" << std::endl;
    retval |= 1<<5;
  }


  return retval;
}
