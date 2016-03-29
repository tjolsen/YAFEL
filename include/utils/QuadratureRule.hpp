#ifndef _YAFEL_QUADRATURERULE_HPP
#define _YAFEL_QUADRATURERULE_HPP

/*
 * Class defining the interface for quadrature rules.
 * Uses CRTP for static polymorphism. Users will instantiate
 * children of this class that implement specific quadrature rules.
 *
 * Children of this class will include TensorProductQuadrature, (...others for simplices)
 *
 */


#include "yafel_globals.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class QuadratureRule {

public:
  using coordinate_type = Tensor<NSD,1,double>;
  using size_type = typename Tensor<NSD,1,double>::size_type;
  using value_type = typename Tensor<NSD,1,double>::value_type;

  QuadratureRule(size_type nPoints) : nodes(nPoints), weights(nPoints) {}
  
  virtual value_type weight(size_type qpi) const {return weights[qpi];}
  virtual coordinate_type qp(size_type qpi) const {return nodes[qpi];}
  inline size_type n_qp() const {return weights.size();}

  std::vector<coordinate_type> nodes;
  std::vector<double> weights;
};

YAFEL_NAMESPACE_CLOSE

#endif
