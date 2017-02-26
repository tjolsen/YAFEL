//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_SEQUENCE_FUNCTIONS_HPP
#define YAFEL_SEQUENCE_FUNCTIONS_HPP

#include "sequences.hpp"
#include <stdexcept>

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

// get the "idx"-th value in a sequence.
constexpr auto index_at(sequence<>, int idx) {
    // goofy error handling to make code choke at compile-time if "idx" was a bad value (ie, too big/negative)
    return (idx==0 || idx) ? throw std::invalid_argument("Bad value in index_at") : -1;
}

template<int S, int ...SS>
constexpr auto index_at(sequence<S, SS...>, int idx) {
    return (idx==0) ? S : index_at(sequence<SS...>(), idx-1);
}

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_SEQUENCE_FUNCTIONS_HPP
