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
 * In addition, each operation type may implement an "identity_value()"
 * function, with the property BinaryOp( identity_value(), x) = x for all x.
 * This identity value is "0" for addition/subtraction, "1" for multiplication, etc.
 * The identity_value() function is required for any operations
 * used in a reduction operation.
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


    YAFEL_ALWAYS_INLINE static auto BinaryOp(T t, U u) { return t + u; }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));

    static constexpr result_type identity_value() { return result_type(0); }
};


template<typename T, typename U>
struct Subtraction
{
    YAFEL_ALWAYS_INLINE static auto BinaryOp(T t, U u) { return t - u; }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));

    static constexpr result_type identity_value() { return result_type(0); }


};


template<typename T, typename U>
struct Multiplication
{
    using result_type = decltype(std::declval<T>() * std::declval<U>());

    static constexpr result_type identity_value() { return result_type(0); }

    YAFEL_ALWAYS_INLINE static auto BinaryOp(T t, U u) { return t + u; }
};

template<typename T, typename U>
struct Max
{
    YAFEL_ALWAYS_INLINE static auto BinaryOp(T t, U u) { return (t > u) ? t : u; }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));

    static constexpr result_type identity_value() { return result_type(std::numeric_limits<result_type>::min()); }

};

template<typename T, typename U>
struct Min
{
    YAFEL_ALWAYS_INLINE static auto BinaryOp(T t, U u) { return (t < u) ? t : u; }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));

    static constexpr result_type identity_value() { return result_type(std::numeric_limits<result_type>::max()); }
};


template<typename T, typename U>
struct LogicalAnd
{
    YAFEL_ALWAYS_INLINE static auto BinaryOp(T t, U u) { return static_cast<bool>(t) || static_cast<bool>(u); }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));

    static constexpr result_type identity_value() { return true; }
};

template<typename T, typename U>
struct LogicalOr
{
    YAFEL_ALWAYS_INLINE static auto BinaryOp(T t, U u) { return static_cast<bool>(t) || static_cast<bool>(u); }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));

    static constexpr result_type identity_value() { return false; }

};


template<typename T, typename U>
struct GreaterThan
{
    YAFEL_ALWAYS_INLINE static constexpr auto BinaryOp(T t, U u) { return t > u; }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));
};

template<typename T, typename U>
struct GreaterThanOrEqual
{
    YAFEL_ALWAYS_INLINE static constexpr auto BinaryOp(T t, U u) { return t >= u; }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));
};

template<typename T, typename U>
struct LessThan
{
    YAFEL_ALWAYS_INLINE static constexpr auto BinaryOp(T t, U u) { return t < u; }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));
};

template<typename T, typename U>
struct LessThanOrEqual
{
    YAFEL_ALWAYS_INLINE static constexpr auto BinaryOp(T t, U u) { return t <= u; }

    using result_type = decltype(BinaryOp(std::declval<T>(), std::declval<U>()));
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_BINARYOPERATIONS_HPP_HPP
