#ifndef __YAFEL_SPATIALFUNCTION_HPP
#define __YAFEL_SPATIALFUNCTION_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/Tensor.hpp"

#include <functional>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD, typename R, typename dataType=double>
using SpatialFunction = std::function<R(const Tensor<NSD,1,dataType>&)>;


YAFEL_NAMESPACE_CLOSE

#endif
