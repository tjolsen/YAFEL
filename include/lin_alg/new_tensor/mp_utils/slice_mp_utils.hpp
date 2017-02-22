//
// Created by tyler on 2/22/17.
//

#ifndef YAFEL_SLICE_MP_UTILS_HPP
#define YAFEL_SLICE_MP_UTILS_HPP

#include "yafel_globals.hpp"
#include "sequences.hpp"
#include "sequence_functions.hpp"
#include "TypeList.hpp"

YAFEL_NAMESPACE_OPEN

/**
 * \class slice_sentinel
 * \brief Special struct to be used to indicate a TensorSlice operation
 */
struct slice_sentinel
{
};



//-----------------------------------------------------------------
//make parent slice sequence
//-----------------------------------------------------------------

template<int S>
constexpr auto slice_stride_helper(std::true_type, sequence<S>)
{
    return sequence<S>();
}

template<int S>
constexpr auto slice_stride_helper(std::false_type, sequence<S>)
{
    return sequence<>();
}


template<typename Arg, int S>
constexpr auto make_slice_strides(sequence<S>, Arg)
{
    return slice_stride_helper(typename std::is_same<Arg, slice_sentinel>::type(), sequence<S>());
};

template<typename Arg, typename ...Args,
        int S, int ...SS>
constexpr auto make_slice_strides(sequence<S, SS...>, Arg, Args ...args)
{
    return seq_cat(slice_stride_helper(typename std::is_same<Arg, slice_sentinel>::type(), sequence<S>()),
                   make_slice_strides(sequence<SS...>(), args...));
};


//------------------------------------------------------------------------
// Compute slice offset
//------------------------------------------------------------------------
template<typename INT, int S>
constexpr int slice_offset_helper(std::true_type, sequence<S>, INT)
{
    return 0;
};


template<typename INT, int S, typename=typename std::enable_if<std::is_integral<INT>::value>::type>
constexpr int slice_offset_helper(std::false_type, sequence<S>, INT idx)
{
    return S*idx;
}


template<typename Arg, int S>
constexpr int make_slice_offset(sequence<S>, Arg arg)
{
    return slice_offset_helper(typename std::is_same<Arg,slice_sentinel>::type(),
                               sequence<S>(), arg);
};



template<typename Arg, typename ...Args,
        int S, int ...SS>
constexpr int make_slice_offset(sequence<S,SS...>, Arg arg, Args ...args)
{
    return slice_offset_helper(typename std::is_same<Arg,slice_sentinel>::type(),
                               sequence<S>(), arg)
           + make_slice_offset(sequence<SS...>(), args...);
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_SLICE_MP_UTILS_HPP
