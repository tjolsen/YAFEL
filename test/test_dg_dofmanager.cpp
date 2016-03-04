#include "yafel_globals.hpp"
#include "mesh/RectilinearMesh.hpp"
#include "utils/DG_DoFManager.hpp"

using namespace yafel;


bool test_1() {

  RectilinearMesh<2> M;
  DG_DoFManager<RectilinearMesh<2>,2> dofm(M,1);
  
  bool good = true;
  
  for(std::size_t i=0; i<M.n_elements(); ++i) {
    for(std::size_t n=0; n<M.element(i).size(); ++n) {
      
      good = good && (i*4 + n) == dofm.global_index(i,n,0);
      
    }
  }
  
  return good;
}

bool test_2() {

  RectilinearMesh<2> M;
  DG_DoFManager<RectilinearMesh<2>,2> dofm(M,2);
  
  bool good = true;
  
  for(std::size_t i=0; i<M.n_elements(); ++i) {
    for(std::size_t n=0; n<M.element(i).size(); ++n) {
      for(std::size_t j=0; j<2; ++j) {
        good = good && (i*4*2 + n*2 + j) == dofm.global_index(i,n,j);
      }
    }
  }
  
  return good;
}


int main() {

  int retval(0);

  if(!test_1()) {
    retval |= 1<<0;
    std::cerr << "Failed test_1()" << std::endl;
  }
  if(!test_2()) {
    retval |= 1<<1;
    std::cerr << "Failed test_2()" << std::endl;
  }


  return retval;
}
