#include "yafel.hpp"
#include <iostream>
#include <cstdlib>

using namespace yafel;


int main(int argc, char **argv) {
  
  unsigned npts = 3;
  if(argc >=2)
    npts = atoi(argv[1]);
  
  
  GaussLobattoQuadrature GQ(npts);
  double wsum = 0;
  for(unsigned i=0; i<GQ.get_nqp(); ++i) {
    printf("x = %f\tw = %f\n", GQ.node(i), GQ.weight(i));
    wsum += GQ.weight(i);
  }
  
  printf("sum(w) = %f\n", wsum);

  return 0;
}
