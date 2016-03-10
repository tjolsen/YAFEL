#include "utils/GaussLegendreQuadrature.hpp"

YAFEL_NAMESPACE_OPEN


template<>
GaussLegendreQuadrature<1>::GaussLegendreQuadrature(size_type Npoints) :
QuadratureRule<1>(Npoints)
{
  std::vector<value_type> n1d(Npoints, 0);
  std::vector<value_type> w1d(Npoints, 0);
  compute_1d_points(n1d, w1d, Npoints);

  for(size_type i=0; i<Npoints; ++i) {
    this->weights[i] = w1d[i];
    this->nodes[i] = coordinate_type{n1d[i]};
  }
}

template<>
GaussLegendreQuadrature<2>::GaussLegendreQuadrature(size_type Npoints) :
QuadratureRule<2>(Npoints*Npoints)
{
  std::vector<value_type> n1d(Npoints, 0);
  std::vector<value_type> w1d(Npoints, 0);
  compute_1d_points(n1d, w1d, Npoints);

  for(size_type i=0; i<Npoints; ++i) {
    for(size_type j=0; j<Npoints; ++j) {
      this->weights[i*Npoints + j] = w1d[i]*w1d[j];
      this->nodes[i*Npoints + j] = coordinate_type{n1d[j], n1d[i]};
    }
  }
}

template<>
GaussLegendreQuadrature<3>::GaussLegendreQuadrature(size_type Npoints) :
QuadratureRule<3>(Npoints*Npoints*Npoints)
{
  std::vector<value_type> n1d(Npoints, 0);
  std::vector<value_type> w1d(Npoints, 0);
  compute_1d_points(n1d, w1d, Npoints);

  for(size_type i=0; i<Npoints; ++i) {
    for(size_type j=0; j<Npoints; ++j) {
      for(size_type k=0; k<Npoints; ++k) {
        this->weights[(i*Npoints + j)*Npoints + k] = w1d[i]*w1d[j]*w1d[k];
        this->nodes[(i*Npoints + j)*Npoints + k] = coordinate_type{n1d[k], n1d[j], n1d[i]};
      }
    }
  }
}


YAFEL_NAMESPACE_CLOSE
