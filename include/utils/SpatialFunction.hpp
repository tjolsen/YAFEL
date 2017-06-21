#ifndef _YAFEL_SPATIALFUNCTION_HPP
#define _YAFEL_SPATIALFUNCTION_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/tensors.hpp"

#include <functional>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD, typename R, typename dataType=double>
using SpatialFunction = std::function<R(const Tensor<NSD,1,dataType>&)>;


YAFEL_NAMESPACE_CLOSE

#endif
