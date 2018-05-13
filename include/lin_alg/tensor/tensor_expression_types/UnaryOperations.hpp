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
    static constexpr auto UnaryOp(T t) { return -t; }

    using result_type = decltype(UnaryOp(std::declval<T>()));
};


template<typename T>
struct Sqrt
{
    static auto UnaryOp(T t)
    {
        using std::sqrt;
        return sqrt(t);
    }

    using result_type = decltype(UnaryOp(std::declval<T>()));

};


template<typename T>
struct Sign
{
    static auto UnaryOp(T t) { return (t >= 0) ? 1 : -1; }

    using result_type = decltype(UnaryOp(std::declval<T>()));
};


template<typename T>
struct Sin
{
    static auto UnaryOp(T t)
    {
        using std::sin;
        return sin(t);
    }

    using result_type = decltype(UnaryOp(std::declval<T>()));
};

template<typename T>
struct Cos
{
    static auto UnaryOp(T t)
    {
        using std::cos;
        return cos(t);
    }

    using result_type = decltype(UnaryOp(std::declval<T>()));
};

template<typename T>
struct Round
{
    static auto UnaryOp(T t)
    {
        using std::round;
        return round(t);
    }

    using result_type = decltype(UnaryOp(std::declval<T>()));
};


template<typename T>
struct Floor {
    static auto UnaryOp(T t) {
        using std::floor;
        return floor(t);
    }
    using result_type = decltype(UnaryOp(std::declval<T>()));
};

template<typename T>
struct Ceil {
    static auto UnaryOp(T t) {
        using std::ceil;
        return ceil(t);
    }
    using result_type = decltype(UnaryOp(std::declval<T>()));
};

template<typename T>
struct Abs
{
    static auto UnaryOp(T t)
    {
        using std::abs;
        return abs(t);
    }

    using result_type = decltype(UnaryOp(std::declval<T>())); //decltype(abs(std::declval<T>()));
};

template<typename T>
struct IsNan
{
    using result_type = bool;

    static auto UnaryOp(T t)
    {
        using std::isnan;
        return isnan(t);
    }
};

template<typename T>
struct TensorCaster {
    template<typename U>
    struct Cast {
        static auto UnaryOp(U u) {

        }
        using result_type = T;
    };
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_UNARYOPERATIONS_HPP
