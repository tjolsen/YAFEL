#include "utils/GaussLobattoQuadrature.hpp"


YAFEL_NAMESPACE_OPEN

/*
 * Specialized Ctors for each NSD to prevent need for recusion
 */
template<>
GaussLobattoQuadrature<1>::GaussLobattoQuadrature(size_type Npoints)
: QuadratureRule<1>(Npoints), nodes_1d(Npoints), pairs(Npoints)
{
  std::vector<value_type> w1d(Npoints,0);
  compute_1d_points(nodes_1d, w1d, Npoints);

  for(size_type i=0; i<Npoints; ++i) {
    this->weights[i] = w1d[i];
    this->nodes[i] = coordinate_type{nodes_1d[i]};
    this->pairs[i] = Tensor<1,1,size_type>{i};
  }
}

template<>
GaussLobattoQuadrature<2>::GaussLobattoQuadrature(size_type Npoints)
: QuadratureRule<2>(Npoints*Npoints), nodes_1d(Npoints), pairs(Npoints*Npoints)
{
  std::vector<value_type> w1d(nodes.size(),0);
  compute_1d_points(nodes_1d, w1d, Npoints);

  for(size_type i=0; i<Npoints; ++i) {
    for(size_type j=0; j<Npoints; ++j) {
      this->weights[i*Npoints + j] = w1d[i]*w1d[j];
      this->nodes[i*Npoints + j] = coordinate_type{nodes_1d[j], nodes_1d[i]};
      this->pairs[i*Npoints + j] = Tensor<2,1,size_type>{j,i};
    }
  }
}

template<>
GaussLobattoQuadrature<3>::GaussLobattoQuadrature(size_type Npoints)
: QuadratureRule<3>(Npoints*Npoints*Npoints), nodes_1d(Npoints), pairs(Npoints*Npoints*Npoints)
{
  std::vector<value_type> w1d(nodes.size(),0);
  compute_1d_points(nodes_1d, w1d, Npoints);

  for(size_type i=0; i<Npoints; ++i) {
    for(size_type j=0; j<Npoints; ++j) {
      for(size_type k=0; k<Npoints; ++k) {
        this->weights[(i*Npoints + j)*Npoints + k] = w1d[i]*w1d[j]*w1d[k];
        this->nodes[(i*Npoints + j)*Npoints + k] = coordinate_type{nodes_1d[k], nodes_1d[j], nodes_1d[i]};
        this->pairs[(i*Npoints + j)*Npoints + k] = Tensor<3,1,size_type>{k,j,i};
      }
    }
  }
}


YAFEL_NAMESPACE_CLOSE
