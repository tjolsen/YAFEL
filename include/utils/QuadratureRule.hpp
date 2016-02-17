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


#include "yafel_gloabls.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

template<typename T>
class QuadratureRule {

public:
  using coordinate_type = Tensor<NSD,1,double>;
  using size_type = typename Tensor<NSD,1,double>::size_type;
  using value_type = typename Tensor<NSD,1,double>::value_type;
  
  value_type weight(size_type qpi) const {return static_cast<T const&>(*this).weight(qpi);}
  coordinate_type qp(size_type qpi) const {return static_cast<T const&>(*this).qp(qpi);}
  size_type n_qp() const {return static_cast<T const&>(*this).n_qp();}


  operator T&() {return static_cast<T&>(*this);}
  operator T const&() const {return static_cast<T const&>(*this);}
  
};

YAFEL_NAMESPACE_CLOSE

#endif
