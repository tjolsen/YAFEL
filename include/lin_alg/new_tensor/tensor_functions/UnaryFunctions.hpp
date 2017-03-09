//
// Created by tyler on 3/9/17.
//

#ifndef YAFEL_UNARYFUNCTIONS_HPP
#define YAFEL_UNARYFUNCTIONS_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorCwiseUnaryOp.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/UnaryOperations.hpp"

YAFEL_NAMESPACE_OPEN

template<typename T, int D, int R, typename dt, bool b>
auto sqrt(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,Sqrt>(te);
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_UNARYFUNCTIONS_HPP
