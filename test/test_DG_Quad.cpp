#include "yafel_globals.hpp"
#include "element/DG_Quad.hpp"
#include "mesh/RectilinearMesh.hpp"
#include "utils/GaussLegendreQuadrature.hpp"
#include "utils/GaussLobattoQuadrature.hpp"
#include "utils/DG_DoFManager.hpp"
#include "utils/DualNumber.hpp"
#include "lin_alg/tensor/Tensor.hpp"

#include <iostream>

using namespace yafel;


// create 1-element mesh on [0, 1]^2 and get the first detJ. (should equal 1/4)
// tons of stuff has to go right (and has already been debugged en route to this test)
// in order to get this correct.
bool test_1() {
  
  RectilinearMesh<2> M(std::vector<double>{1,1}, std::vector<std::size_t>{1,1});

  //use polynomials of order 15 along each edge, and 15-point 1d integration rule (super overkill!)
  //hopefully this high-order test exposes numerical stability issues, if present (which they are not)
  std::size_t N = 15;
  GaussLegendreQuadrature<2> Q2(N);
  GaussLegendreQuadrature<1> Q1(N);
  DG_DoFManager<RectilinearMesh<2>,2> dofm(M, 1);

  DG_Quad<> DGQ(N, dofm, Q2, Q1);

  DGQ.update_element(M, 0);
  
  return std::abs(DGQ.detJ[0] - 1.0/4.0) < 1.0e-8;
}


/*
 * use shape_gradients to compute the gradient of a function (x^2 + y^2)
 * and compare with result obtained by using automatic differentiation
 */
template<typename T>
T func(T x, T y) {
  return x*x + y*y;
}

bool test_2() {

  RectilinearMesh<2> M(std::vector<double>{1,1}, std::vector<std::size_t>{1,1});

  std::size_t N = 2;
  GaussLegendreQuadrature<2> Q2(N);
  GaussLegendreQuadrature<1> Q1(N);
  DG_DoFManager<RectilinearMesh<2>,2> dofm(M, 1);

  DG_Quad<> DGQ(N, dofm, Q2, Q1);

  DGQ.update_element(M, 0);
  
  bool good = true;
  for(std::size_t qpi=0; qpi<DGQ.Q2D.n_qp(); ++qpi) {
    
    Tensor<2,1,double> xqp = DGQ.xval(DGQ.Q2D.qp(qpi));

    double f_x = func(DualNumber<double>(xqp(0),1), DualNumber<double>(xqp(1),0)).second;
    double f_y = func(DualNumber<double>(xqp(0),0), DualNumber<double>(xqp(1),1)).second;

    Tensor<2,1,double> gdual{f_x, f_y};
    Tensor<2,1,double> g;    
    for(std::size_t A=0; A<DGQ.nodes_per_element; ++A) {
      auto x = DGQ.nodes_x[A];
      g(0) += func(x(0), x(1))*DGQ.shape_gradients[qpi](A,0);
      g(1) += func(x(0), x(1))*DGQ.shape_gradients[qpi](A,1);
       
   }

    double normdiff = contract<1>(g-gdual,g-gdual);
    good = good &&  normdiff < 1.0e-8;
  }
  
  
  /*
    for(std::size_t qpi=0; qpi<DGQ.Q2D.n_qp(); ++qpi) {
    
    for(std::size_t A=0; A<DGQ.nodes_per_element; ++A) {
    for(std::size_t i=0; i<2; ++i) {
    std::cout << DGQ.shape_gradients[qpi](A,i) << "\t";
    }
    std::cout << std::endl;
    }
    std::cout << std::endl;
    }
  */
  
  return good;
}

int main() {

  int retval = 0;

  if(!test_1()) {
    std::cout << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }
  if(!test_2()) {
    std::cout << "Failed test_2()" << std::endl;
    retval |= 1<<1;
  }

  return retval;
}
