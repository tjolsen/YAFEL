//
// Created by tyler on 3/8/17.
//

#ifndef YAFEL_BINARYOPERATIONS_HPP_HPP
#define YAFEL_BINARYOPERATIONS_HPP_HPP

/**
 * \file
 *
 * This file defines a set of binary operations that are
 * passed in as template arguments to TensorCwiseBinaryOp.
 * They implment the binary operation associated with their
 * name, and are meant to cut down on the redundant definitions
 * of many different TensorExpressions for component-wise binary
 * operations.
 *
 * Some operations are "builtin" functions: (+, -, *[scalar], /[scalar])
 * Other operations may be functions that operate component-wise
 * (max, min, pow, ...)
 *
 * Each operation type must implement a "BinaryOp(T t, U u)" function,
 * and the rest is taken care of by TensorCwiseBinaryOp.
 *
 * A similar file exists for unary operations.
 */


#include "yafel_globals.hpp"
#include "utils/ScalarTraits.hpp"
#include <type_traits>

YAFEL_NAMESPACE_OPEN

template<typename T, typename U>
struct Addition
{
    using result_type = decltype(T() + U());

    static auto BinaryOp(T t, U u) { return t + u; }
};


template<typename T, typename U>
struct Subtraction
{
    using result_type = decltype(T() - U());

    static auto BinaryOp(T t, U u) { return t - u; }
};


template<typename T, typename U>
struct Max
{
    using result_type = decltype(T() * U());

    static auto BinaryOp(T t, U u) { return (t > u) ? t : u; }
};

template<typename T, typename U>
struct Min
{
    using result_type = decltype(T() * U());

    static auto BinaryOp(T t, U u) { return (t < u) ? t : u; }
};



YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_BINARYOPERATIONS_HPP_HPP
