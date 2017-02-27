//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_SEQUENCES_HPP
#define YAFEL_SEQUENCES_HPP

#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

template<int ...S>
struct sequence{
    constexpr static int size() {return sizeof...(S);}
};

//----------------------------------------------------------
//Generate a sequence of N integer values "V"
template<int V, int N, int ...SS>
struct value_sequence : value_sequence<V,N-1,V, SS...> {};

template<int V, int ...SS>
struct value_sequence<V,0,SS...>
{
    using type = sequence<SS...>;
};

//----------------------------------------------------------
template<int N, int ...S>
struct counting_sequence : public counting_sequence<N-1, N-1, S...>{};

template<int ...S>
struct counting_sequence<0,S...> {
    using type = sequence<S...>;
};



//----------------------------------------------------------
template<int N, int Base, int S, int ...SS>
struct geometric_sequence : public geometric_sequence<N-1, Base, Base*S, S, SS...>{};

template<int Base, int S, int ...SS>
struct geometric_sequence<0,Base,S,SS...>
{
    using type = sequence<SS...>;
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_SEQUENCES_HPP
