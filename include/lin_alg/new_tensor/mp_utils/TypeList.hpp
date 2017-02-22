//
// Created by tyler on 2/22/17.
//

#ifndef YAFEL_TYPELIST_HPP
#define YAFEL_TYPELIST_HPP

#include "yafel_globals.hpp"
#include <type_traits>



YAFEL_NAMESPACE_OPEN


/**
 * \brief Template container to hold a list of types
 * @tparam Args
 */
template<typename ...Args>
struct type_list {};


// Test for all type T in type list
template<typename T>
constexpr bool all_of_type(type_list<>) {
    return true;
}

template<typename T, typename Arg, typename ...Args>
constexpr bool all_of_type(type_list<Arg, Args...>)
{
    return std::is_same<T,Arg>::value && all_of_type<T>(type_list<Args...>());
};



// Test if type T contained in type list
template<typename T>
constexpr bool contains(type_list<>) {
    return false;
};

template<typename T, typename Arg, typename ...Args>
constexpr bool contains(type_list<Arg, Args...>) {
    return std::is_same<T,Arg>::value || contains<T>(type_list<Args...>());
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TYPELIST_HPP
