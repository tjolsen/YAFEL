//
// Created by tyler on 3/9/17.
//

#ifndef YAFEL_BINARYFUNCTIONS_HPP
#define YAFEL_BINARYFUNCTIONS_HPP

/**
 * \file
 * \brief Define component-wise binary functions on tensors
 */

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/BinaryOperations.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorCwiseBinaryOp.hpp"


YAFEL_NAMESPACE_OPEN

/**
 * \brief Compute component-by-component max
 */
template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b1, bool b2>
auto max(const TensorExpression<T1,D,R,dt1,b1> &lhs,
         const TensorExpression<T1,D,R,dt2,b2> &rhs)
{
    return TensorCwiseBinaryOp<T1,T2,D,R,dt1,dt2,b1,b2,Max>(lhs,rhs);
}

/**
 * \brief Compute component-by-component min
 */
template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b1, bool b2>
auto min(const TensorExpression<T1,D,R,dt1,b1> &lhs,
         const TensorExpression<T1,D,R,dt2,b2> &rhs)
{
    return TensorCwiseBinaryOp<T1,T2,D,R,dt1,dt2,b1,b2,Min>(lhs,rhs);
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_BINARYFUNCTIONS_HPP
