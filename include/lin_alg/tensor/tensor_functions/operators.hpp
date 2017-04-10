//
// Created by tyler on 3/8/17.
//

#ifndef YAFEL_OPERATORS_HPP
#define YAFEL_OPERATORS_HPP

/**
 * \file
 *
 * \brief This file defines operators that operate on TensorExpressions
 */

#include "yafel_globals.hpp"
#include "lin_alg/tensor/tensor_expression_types/BinaryOperations.hpp"
#include "lin_alg/tensor/tensor_expression_types/UnaryOperations.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorCwiseBinaryOp.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorCwiseUnaryOp.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorScaled.hpp"

YAFEL_NAMESPACE_OPEN

//-------------------------------------------------------------------
//              Component-wise Binary Operations
//-------------------------------------------------------------------


// Addition
template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b1, bool b2>
auto operator+(const TensorExpression<T1, D, R, dt1, b1> &lhs,
               const TensorExpression<T2, D, R, dt2, b2> &rhs)
{
    return TensorCwiseBinaryOp<T1, T2, D, R, dt1, dt2, b1, b2, Addition>(lhs, rhs);
}


// Subtraction
template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b1, bool b2>
auto operator-(const TensorExpression<T1, D, R, dt1, b1> &lhs,
               const TensorExpression<T2, D, R, dt2, b2> &rhs)
{
    return TensorCwiseBinaryOp<T1, T2, D, R, dt1, dt2, b1, b2, Subtraction>(lhs, rhs);
}


// Multiplication by scalar (from rhs)
template<typename T1, int D, int R, typename dt, bool b, typename U,
        typename=typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator*(const TensorExpression<T1, D, R, dt, b> &lhs, const U &rhs)
{
    return TensorScaled<T1, U, D, R, dt>(lhs, rhs);
}


// Multiplication by scalar (from lhs)
template<typename T1, int D, int R, typename dt, bool b, typename U,
        typename=typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator*(const U &lhs, const TensorExpression<T1, D, R, dt, b> &rhs)
{
    return TensorScaled<T1, U, D, R, dt>(rhs, lhs);
}

// Division (recast as multiplication by inverse, which is faster)
template<typename T1, int D, int R, typename dt, bool b, typename U,
        typename=typename std::enable_if<std::is_arithmetic<U>::value>::type>
auto operator/(const TensorExpression<T1, D, R, dt, b> &lhs, const U &rhs)
{
    using result_type = decltype(dt() / U());
    return TensorScaled<T1, result_type, D, R, dt>(lhs, result_type(1) / rhs);
}

//-------------------------------------------------------------------
//              Component-wise Unary Operations
//-------------------------------------------------------------------
template<typename T1, int D, int R, typename dt1, bool b1>
auto operator-(const TensorExpression<T1, D, R, dt1, b1> &te)
{
    return TensorCwiseUnaryOp<T1, D, R, dt1, b1, Negation>(te);
}


YAFEL_NAMESPACE_CLOSE
#endif //YAFEL_OPERATORS_HPP
