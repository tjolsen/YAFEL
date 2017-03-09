//
// Created by tyler on 3/9/17.
//

#ifndef YAFEL_SCALARTRAITS_HPP
#define YAFEL_SCALARTRAITS_HPP

#include "yafel_globals.hpp"
#include <type_traits>

YAFEL_NAMESPACE_OPEN

template<typename T>
struct ScalarTraits
{
    static constexpr bool isYafelScalar() { return std::is_arithmetic<T>::value; }
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_SCALARTRAITS_HPP
