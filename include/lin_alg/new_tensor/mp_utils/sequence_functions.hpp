//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_SEQUENCE_FUNCTIONS_HPP
#define YAFEL_SEQUENCE_FUNCTIONS_HPP

#include "sequences.hpp"
#include <type_traits>

YAFEL_NAMESPACE_OPEN


// Sequence sum
constexpr int sum(sequence<>) { return 0; }

template<int S, int ...SS>
constexpr int sum(sequence<S, SS...>) { return S + sum(sequence<SS...>()); }


// Concatenate sequences
template<int ...PP, int ...SS>
constexpr auto seq_cat(sequence<PP...>, sequence<SS...>)
{
    return sequence<PP..., SS...>();
}


//reverse sequences
constexpr auto seq_reverse(sequence<>)
{
    return sequence<>();
}

template<int S, int ...SS>
constexpr auto seq_reverse(sequence<S,SS...>)
{
    return seq_cat(seq_reverse(sequence<SS...>()), sequence<S>());
};

// get the "idx"-th value in a sequence.
template<int I, int S, int ...SS, typename=typename std::enable_if<0 < sizeof...(SS)>::type>
constexpr auto index_at(sequence<S, SS...>, sequence<I>)
{
    //static_assert(I<=sizeof...(SS),"index_at: 'I' too large");
    return (I == 0) ? S : index_at(sequence<SS...>(), sequence<I - 1>());
}

template<int I, int S>
constexpr auto index_at(sequence<S>, sequence<I>)
{
    // Must allow negative 'I' in order to instantiate function chain
    static_assert(I <= 0, "index_at: base case error. 'I' likely too large.");
    return S;
};


// Test if all elements are >= some other value
template<int I>
constexpr bool all_ge(sequence<>, sequence<I>) { return true; }

template<int I, int S, int ...SS>
constexpr bool all_ge(sequence<S, SS...>, sequence<I>)
{
    return (S >= I) && all_ge(sequence<SS...>(), sequence<I>());
};


// Test if all elements are < some other value
template<int I>
constexpr bool all_lt(sequence<>, sequence<I>) { return true; }

template<int I, int S, int ...SS>
constexpr bool all_lt(sequence<S, SS...>, sequence<I>)
{
    return (S < I) && all_lt(sequence<SS...>(), sequence<I>());
};


//Test if sequence contains an element
template<int I>
constexpr bool seq_contains(sequence<>, sequence<I>) { return false; }

template<int I, int S, int ...SS>
constexpr bool seq_contains(sequence<S, SS...>, sequence<I>)
{
    return I == S || seq_contains(sequence<SS...>(), sequence<I>());
};

//Test if sequence contains a duplicate (quadratic version)
template<int S>
constexpr bool contains_duplicates(sequence<S>) { return false; }

template<int S, int ...SS, typename=typename std::enable_if<0 < sizeof...(SS)>::type>
constexpr bool contains_duplicates(sequence<S, SS...>)
{
    return seq_contains(sequence<SS...>(), sequence<S>())
           || contains_duplicates(sequence<SS...>());
}


//Find element in sequence (first occurrence). Return sequence::size() if not in sequence.
template<int Val>
constexpr int seq_find(sequence<Val>, sequence<>) { return 0; }

template<int Val, int S, int ...SS>
constexpr int seq_find(sequence<Val>, sequence<S, SS...>)
{
    return (Val == S) ? 0 : 1 + seq_find(sequence<Val>(), sequence<SS...>());
};

// Get last N elements of sequence
template<int N, int remaining, int S, int ...SS>
struct last_N : public last_N<N, remaining - 1, SS...>
{
};

template<int N, int ...SS>
struct last_N<N, N, SS...>
{
    using type = sequence<SS...>;
};

template<int ...SS>
constexpr auto get_last_N(sequence<0>,sequence<SS...>)
{
    return sequence<>();
}

template<int N, int ...SS>
constexpr auto get_last_N(sequence<N>, sequence<SS...>)
{
    static_assert(N<=sizeof...(SS),"Error: Cannot get last N of sequence. Sequence too short.");
    return typename last_N<N, sizeof...(SS), SS...>::type();
};


// Get first N elements of sequence
template<int N, int ...SS>
constexpr auto get_first_N(sequence<N>, sequence<SS...>)
{
    static_assert(N <= sizeof...(SS),"Error: cannot get first N of sequence. Sequence too short");

    return seq_reverse(get_last_N(sequence<N>(), seq_reverse(sequence<SS...>())));
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_SEQUENCE_FUNCTIONS_HPP
