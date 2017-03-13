//
// Created by tyler on 3/12/17.
//

#ifndef YAFEL_MESH_TYPEDEFS_HPP
#define YAFEL_MESH_TYPEDEFS_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/Tensor.hpp"

YAFEL_NAMESPACE_OPEN

//Unified coordinate<T> type across library.
// No longer templated on the "NSD" parameter.
template<typename T=double>
using coordinate = Tensor<3, 1, T>;

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_MESH_TYPEDEFS_HPP
