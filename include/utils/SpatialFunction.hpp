#ifndef __YAFEL_SPATIALFUNCTION_HPP
#define __YAFEL_SPATIALFUNCTION_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/Tensor.hpp"

YAFEL_NAMESPACE_OPEN

template<unsigned NSD, typename T=double>
class SpatialFunction {

private:
  T (*func)(const Tensor<NSD,1,T> &x);

public:
  SpatialFunction(T (*fp)(const Tensor<NSD,1,T> &x)) : func(fp) {}
  T operator()(const Tensor<NSD,1,T> &x) const {return func(x);}
};

YAFEL_NAMESPACE_CLOSE

#endif
