#include "yafel_globals.hpp"
#include "element/DG_Quad.hpp"
#include "mesh/RectilinearMesh.hpp"
#include "utils/GaussLegendreQuadrature.hpp"
#include "utils/GaussLobattoQuadrature.hpp"
#include "utils/DG_DoFManager.hpp"
#include "lin_alg/tensor/Tensor.hpp"

#include <iostream>

using namespace yafel;


// create 1-element mesh on [0, 1]^2 and get the first detJ. (should equal 1/4)
// tons of stuff has to go right (and has already been debugged en route to this test)
// in order to get this correct.
bool test_1() {
  
  RectilinearMesh<2> M(std::vector<double>{1,1}, std::vector<std::size_t>{1,1});

  //use 15 nodes along each edge, and 15-point 1d integration rule (super overkill!)
  std::size_t N = 15;
  GaussLegendreQuadrature<2> Q2(N);
  GaussLegendreQuadrature<1> Q1(N);
  DG_DoFManager<RectilinearMesh<2>,2> dofm(M, 1);

  DG_Quad<> DGQ(N, dofm, Q2, Q1);

  DGQ.update_element(M, 0);
  
  return std::abs(DGQ.detJ[0] - 1.0/4.0) < 1.0e-8;
}


int main() {

  int retval = 0;

  if(!test_1()) {
    std::cout << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }

  return retval;
}
