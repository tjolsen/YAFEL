#include "yafel_globals.hpp"
#include "mesh/RectilinearMesh.hpp"
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
  
  std::size_t N = 10;
  RectilinearMesh<2> M(std::vector<double>{1,1}, std::vector<std::size_t>{N,N});
  M.build_faces();
  
  //test right neighbor
  for(std::size_t i=0; i<N; ++i) {
    for(std::size_t j=0; j<N; ++j) {
      
      
      
    }
  }
  
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
