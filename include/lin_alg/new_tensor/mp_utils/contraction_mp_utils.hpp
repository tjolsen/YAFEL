//
// Created by tyler on 2/26/17.
//

#ifndef YAFEL_CONTRACTION_MP_UTILS_HPP
#define YAFEL_CONTRACTION_MP_UTILS_HPP

#include "yafel_globals.hpp"
#include "sequences.hpp"
#include "sequence_functions.hpp"

YAFEL_NAMESPACE_OPEN

template<int N, int ...SS>
constexpr auto lhs_parent_strides(sequence<N>, sequence<SS...>)
{
    return get_last_N(sequence<N>(), sequence<SS...>());

}

template<int N, int ...SS>
constexpr auto rhs_parent_strides(sequence<N>, sequence<SS...>)
{
    return get_first_N(sequence<N>(), sequence<SS...>());
}


inline int get_offset(sequence<>, sequence<>, int)
{
    return 0;
}

template<int ...PP>
inline int get_offset(sequence<>, sequence<PP...>, int)
{
    return 0;
}

template<int ...SS>
inline int get_offset(sequence<SS...>, sequence<>, int)
{
    return 0;
}

template<int S, int ...SS, int P, int ...PP>
inline int get_offset(sequence<S, SS...>, sequence<P, PP...>, int idx)
{
    return P * (idx / S) + get_offset(sequence<SS...>(), sequence<PP...>(), idx % S);
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_CONTRACTION_MP_UTILS_HPP
