#include "yafel_globals.hpp"
#include "mesh/RectilinearMesh.hpp"
#include "mesh/GmshMesh.hpp"
#include <vector>
#include <iostream>

using namespace yafel;

bool test_1() {
  
  std::size_t N = 200;
  RectilinearMesh<2> M(std::vector<double>{1,1}, std::vector<std::size_t>{2,2});
  M.build_faces();

  RectilinearMesh<2> MM(std::vector<double>{1,1}, std::vector<std::size_t>{N,N});
  MM.build_faces();

  return M.mesh_faces.size() == 12 && MM.mesh_faces.size()==2*(N+1)*N;
}

bool test_2() {
  
  GmshMesh<2> M("minsquare.msh");
  M.build_faces();
  
  //use Euler's formula for correct number of edges
  //http://nrich.maths.org/1384
  return M.n_elements() - M.mesh_faces.size() + M.n_nodes() == 1;
}


int main() {
  
  int retval(0);

  if(!test_1()) {
    std::cerr << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }
  if(!test_2()) {
    std::cerr << "Failed test_2()" << std::endl;
    retval |= 1<<1;
  }
  
  return retval;
}
