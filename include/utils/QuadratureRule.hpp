#ifndef __YAFEL_QUADRATURERULE_HPP
#define __YAFEL_QUADRATURERULE_HPP

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
  
  virtual value_type weight(size_type qpi) const = 0;
  virtual coordinate_type qp(size_type qpi) const = 0;
  virtual size_type n_qp() const = 0;

};

YAFEL_NAMESPACE_CLOSE

#endif
