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

template<typename T>
struct Negation
{
    using result_type = decltype(-T());

    static constexpr auto UnaryOp(T t) { return -t; }
};


template<typename T>
struct Sqrt
{
    using result_type = decltype(sqrt(std::declval<T>()));

    static auto UnaryOp(T t) { return sqrt(t); }
};


template<typename T>
struct Sign
{
    //always return a signed int from this function.
    using result_type = int;

    static auto UnaryOp(T t) { return (t >= 0) ? 1 : -1; }
};


template<typename T>
struct Sin
{
    using result_type = decltype(sin(std::declval<T>()));

    static auto UnaryOp(T t) { return sin(t); }
};

template<typename T>
struct Cos
{
    using result_type = decltype(sin(std::declval<T>()));

    static auto UnaryOp(T t) { return cos(t); }
};

template<typename T>
struct Round
{
    using result_type = decltype(std::round(std::declval<T>()));

    static auto UnaryOp(T t) { return std::round(t); }
};

template<typename T>
struct Abs
{
    using result_type = decltype(std::abs(std::declval<T>()));
    static auto UnaryOp(T t) { using std::abs; return abs(t); }
};

template<typename T>
struct IsNan
{
    using result_type = bool;
    static auto UnaryOp(T t) { using std::isnan; return isnan(t); }
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_UNARYOPERATIONS_HPP
