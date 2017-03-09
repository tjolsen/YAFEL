//
// Created by tyler on 3/9/17.
//

#ifndef YAFEL_UNARYOPERATIONS_HPP
#define YAFEL_UNARYOPERATIONS_HPP

/**
 * \file
 * \brief This file defines unary component-wise operations on TensorExpressions
 */

#include "yafel_globals.hpp"
#include "utils/ScalarTraits.hpp"
#include <type_traits>
#include <cmath>

YAFEL_NAMESPACE_OPEN

/**
 * \class UnaryOperator
 * \brief Base type of all structs that implement the unary operator policy
 * Can only be instantiated with Yafel Scalar types.
 * This is given by an entry in yafel::ScalarTraits.
 * By default, all builtin arithmetic types (std::is_arithmetic<T>::value == true)
 * are Yafel Scalars.
 * Any other data types (such as DualNumber) must implement their own
 * ScalarTraits class.
 *
 * @tparam T
 */
template<typename T>
        //typename = typename std::enable_if<ScalarTraits<T>::isYafelScalar()>::type>
struct UnaryOperator {};

template<typename T>
struct Negation : UnaryOperator<T>
{
    using result_type = decltype(-T());

    static constexpr auto UnaryOp(T t) { return -t; }
};


template<typename T>
struct Sqrt : UnaryOperator<T>
{
    using result_type = decltype(sqrt(std::declval<T>()));

    static auto UnaryOp(T t) { return sqrt(t); }
};


template<typename T>
struct Sign : UnaryOperator<T>
{
    //always return a signed int from this function.
    using result_type = int;

    static auto UnaryOp(T t) { return (t >= 0) ? 1 : -1; }
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_UNARYOPERATIONS_HPP
