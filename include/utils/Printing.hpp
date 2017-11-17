//
// Created by tyler on 11/16/17.
//

#ifndef YAFEL_PRINTING_HPP
#define YAFEL_PRINTING_HPP

#include "yafel_globals.hpp"
#include "utils/DualNumber.hpp"


#include <iostream>
/**
 * \file
 * This file condenses functions that print yafel data structures
 * into a utility file in order to remove unnecessary <iostream>
 * dependencies throughout the codebase. This should help with
 * some of the slower compile times.
 */

YAFEL_NAMESPACE_OPEN

template<typename T>
std::ostream & operator<<(std::ostream & out, DualNumber<T> x)
{
    out << "{(" << x.first << ") + (" << x.second << ")eps}";
    return out;
}



YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_PRINTING_HPP
