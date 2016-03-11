#include "yafel.hpp"
#include "LinearAdvection.hpp"
#include <vector>
#include <cmath>
#include <string>

using namespace yafel;

int main() {

  Tensor<2,1,double> c0{1,0}; // <-- un-normalized velocity vector. Control direction here.

  std::size_t P = 1;
  GaussLegendreQuadrature<2> Q2D(P+1);
  GaussLegendreQuadrature<1> Q1D(P+1);



  std::vector<double> mesh_dims{2, 1};
  std::vector<std::size_t> dir_elems{50,25};
  Tensor<2,1,double> v_advection = c0/std::sqrt(contract<1>(c0,c0)); // <--normalize to 1 for now
  double dt = .01;
  double T = 1;
  
  SpatialFunction<2,double> u0_func(
                                    [](const Tensor<2,1,double> &x) {
                                      return std::sin(2*3.14159*x(0))*(x(0)<=.5 && x(0)>=0)*(x(1)>=.25 && x(1)<=.75);
                                    }
                                    );

  std::string outfile_base("output");


  AdvectionParameters params(mesh_dims, dir_elems, v_advection, P, Q2D, Q1D, dt, T, u0_func, outfile_base);

  LinearAdvection LA(params);

  LA.run();

  return 0;
}
