#include "yafel_globals.hpp"
#include "output/MatrixVisualization.hpp"
#include "old_handmade_linalg/sparse_coo.hpp"
#include "old_handmade_linalg/sparse_csr.hpp"
#include "old_handmade_linalg/sparse_bcsr.hpp"
#include "mesh/GmshMesh.hpp"

#include <iostream>

using namespace yafel;

// spy an NxN sparse identity matrix. can't check programmatically, must check visually
bool test_1() {
  
  std::size_t N = 100;
  sparse_coo<> A;

  for(std::size_t i=0; i<N; ++i) {
    A.add(i,i,1);
  }
  MatrixVisualization::spy(A);

  return true;
}

// visualize the adjacency of a mesh defined in a gmsh mesh file
// ensure that this works for non-coo matrix type
bool test_2() {
  
  std::string fname = "minsquare.msh";

  GmshMesh<2> M(fname);
  sparse_coo<> A;
  for(std::size_t i=0; i<M.n_elements(); ++i) {

    auto e = M.element(i);
    for(auto a : e) {
      for(auto b : e) {
        A.add(a,b,1);
      }
    }
  }

  sparse_csr<> B(A);
  
  MatrixVisualization::spy(B);

  return true;
}

// give the mesh above a block structure by letting 3 dofs live per node.
// dof_global_num = 3*nodenum + dof_local_num
template<unsigned N>
bool test_3() {
  
  std::string fname = "minsquare.msh";

  GmshMesh<2> M(fname);
  sparse_coo<> A;
  for(std::size_t i=0; i<M.n_elements(); ++i) {

    auto e = M.element(i);
    for(auto a : e) {
      for(auto b : e) {
        
        for(std::size_t dofa=0; dofa<N; ++dofa) {
          for(std::size_t dofb=0; dofb<N; ++dofb) {
            
            A.add(N*a + dofa, N*b + dofb, 1);
            
          }
        }

      }
    }
  }

  sparse_bcsr<N> B(A);
  
  MatrixVisualization::spy(B);

  return true;
}


int main() {

  int retval = 0;
  if(!test_1()) {
    std::cerr << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }
  if(!test_2()) {
    std::cerr << "Failed test_2()" << std::endl;
    retval |= 1<<1;
  }
  if(!test_3<4>()) {
    std::cerr << "Failed test_3<N>()" << std::endl;
    retval |= 1<<2;
  }


  return retval;
}
