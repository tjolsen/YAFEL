//
// Created by tyler on 3/8/17.
//

#ifndef YAFEL_MAP_MP_UTILS_HPP
#define YAFEL_MAP_MP_UTILS_HPP

#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

//---------------------------------------
template<typename T, bool constPtr, bool ptrToConst>
struct map_ptr
{
};

template<typename T>
struct map_ptr<T, true, true>
{
    using type = T const *const;
};

template<typename T>
struct map_ptr<T, true, false>
{
    using type = T *const;
};

template<typename T>
struct map_ptr<T, false, true>
{
    using type = T const *;
};

template<typename T>
struct map_ptr<T, false, false>
{
    using type = T *;
};


//---------------------------------------
template<typename T, bool isConst>
struct map_ref
{
    using type = T &;
};

template<typename T>
struct map_ref<T, true>
{
    using type = T const &;
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_MAP_MP_UTILS_HPP
