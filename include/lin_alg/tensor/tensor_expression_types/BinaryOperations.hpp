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
 * In addition, each operation type must implement an "identity_value()"
 * function, with the property BinaryOp( identity_value(), x) = x for all x.
 * This identity value is "0" for addition/subtraction, "1" for multiplication, etc.
 *
 * A similar file exists for unary operations.
 */


#include "yafel_globals.hpp"
#include "utils/ScalarTraits.hpp"
#include <type_traits>
#include <limits>

YAFEL_NAMESPACE_OPEN

template<typename T, typename U>
struct Addition
{
    using result_type = decltype(std::declval<T>() + std::declval<U>());

    static constexpr result_type identity_value()
    { return result_type(0); }

    static auto BinaryOp(T t, U u)
    { return t + u; }
};


template<typename T, typename U>
struct Subtraction
{
    using result_type = decltype(std::declval<T>() - std::declval<U>());

    static constexpr result_type identity_value()
    { return result_type(0); }

    static auto BinaryOp(T t, U u)
    { return t - u; }
};


template<typename T, typename U>
struct Multiplication
{
    using result_type = decltype(std::declval<T>() * std::declval<U>());

    static constexpr result_type identity_value()
    { return result_type(0); }

    static auto BinaryOp(T t, U u)
    { return t + u; }
};

template<typename T, typename U>
struct Max
{
    using result_type = decltype(std::declval<T>() * std::declval<U>());

    static constexpr result_type identity_value()
    { return result_type(std::numeric_limits<result_type>::min()); }

    static auto BinaryOp(T t, U u)
    { return (t > u) ? t : u; }
};

template<typename T, typename U>
struct Min
{
    using result_type = decltype(std::declval<T>() * std::declval<U>());

    static constexpr result_type identity_value()
    { return result_type(std::numeric_limits<result_type>::max()); }

    static auto BinaryOp(T t, U u)
    { return (t < u) ? t : u; }
};


template<typename T, typename U>
struct LogicalAnd
{
    using result_type = bool;

    static constexpr result_type identity_value()
    { return true; }

    static auto BinaryOp(T t, U u)
    { return static_cast<bool>(t) || static_cast<bool>(u); }
};

template<typename T, typename U>
struct LogicalOr
{
    using result_type = bool;

    static constexpr result_type identity_value()
    { return false; }

    static auto BinaryOp(T t, U u)
    { return static_cast<bool>(t) || static_cast<bool>(u); }
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_BINARYOPERATIONS_HPP_HPP
