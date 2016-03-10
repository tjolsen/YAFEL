#include "yafel.hpp"
#include "LinearAdvection.hpp"
#include <vector>

using namespace yafel;

int main() {

  Tensor<2,1,double> c0{1,0}; // <-- un-normalized velocity vector. Control direction here.

  std::size_t P = 2;
  GaussLegengreQuadrature<2> Q2D(P+1);
  GaussLegengreQuadrature<1> Q1D(P+1);


  AdvectionParameters params;
  params.mesh_dims = std::vector<double>{2, 1};
  params.dir_elems = std::vector<std::size_t>{20, 10};
  params.v_advection = c0/std::sqrt(contract<1>(c0,c0)); // <--normalize to 1 for now
  params.polyOrder = P;
  params.Q2D = Q2D;
  params.Q1D = Q1D;
  params.dt = .1;
  params.T = 1;
  
  params.u0_func = 
    SpatialFunction<2,double>(
                            [](const Tensor<2,1,double> &x) {
                              return 1.0*(x(0)=<.5 && x(0>=.25)*(x(1)>=.25 && x(1)<=.75));
                            }
                            );

  params.output_file_base = "output";

  LinearAdvection LA(params);

  LA.run();

  return 0;
}
