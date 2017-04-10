//
// Created by tyler on 3/14/17.
//

#ifndef YAFEL_YAFEL_TYPEDEFS_HPP
#define YAFEL_YAFEL_TYPEDEFS_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/tensors.hpp"

YAFEL_NAMESPACE_OPEN

//Unified coordinate<T> type across library.
// No longer templated on the "NSD" parameter.
template<typename T=double>
using coordinate = Tensor<3, 1, T>;


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_YAFEL_TYPEDEFS_HPP
