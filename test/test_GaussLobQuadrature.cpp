#include "yafel_globals.hpp"
#include "utils/GaussLobattoQuadrature.hpp"
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace yafel;

template<unsigned NSD>
bool test_1(std::size_t np) {
  GaussLobattoQuadrature<NSD> GQ(np);
  double sum = 0;
  for(auto w : GQ.weights)
    sum += w;
  return std::abs(sum - pow(2.0,NSD)) < 1e-8;
}

int main() {
  
  int retval = 0;
  
  std::size_t npoints = 4;
  if(!test_1<1>(npoints)) {
    std::cerr << "Failed test_1<1>()" << std::endl;
    retval |= 1<<0;
  }
  if(!test_1<2>(npoints)) {
    std::cerr << "Failed test_1<1>()" << std::endl;
    retval |= 1<<0;
  }
  if(!test_1<3>(npoints)) {
    std::cerr << "Failed test_1<1>()" << std::endl;
    retval |= 1<<0;
  }


  return retval;
}
