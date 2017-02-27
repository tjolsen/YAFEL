//
// Created by tyler on 2/26/17.
//

#ifndef YAFEL_MP_TUPLE_UTILS_HPP
#define YAFEL_MP_TUPLE_UTILS_HPP

#include "yafel_globals.hpp"
#include <tuple>

YAFEL_NAMESPACE_OPEN

template<typename T, int N, typename ...TT>
struct homogeneous_tuple : homogeneous_tuple<T,N-1,T,TT...> {};

template<typename T, typename ...TT>
struct homogeneous_tuple<T,0,TT...>
{
    using type = std::tuple<TT...>;
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_MP_TUPLE_UTILS_HPP
