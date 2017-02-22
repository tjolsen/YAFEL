//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_TENSOREXPRESSION_HPP
#define YAFEL_TENSOREXPRESSION_HPP

#include "yafel_globals.hpp"
#include "mp_utils/sequences.hpp"
#include "mp_utils/sequence_functions.hpp"
#include "mp_utils/TypeList.hpp"
#include "mp_utils/slice_mp_utils.hpp"

#include <type_traits>

YAFEL_NAMESPACE_OPEN

// Forward declaration of TensorSlice, so that it can be used from the base class
template<typename TE, int D, int R, typename dataType, bool assignable, int ...PARENT_STRIDES>
class TensorSlice;

template<template<typename, int, int, typename, bool> class TE, typename T, int D, int R, typename dataType, bool assignable, int ...PARENT_STRIDES, typename ...Args>
auto make_slice(const TE<T, D, R, dataType, assignable> &te, Args... args);


/**
 * \class TensorExpression
 * \brief Base class for tensor expressions. All others derive from this.
 * @tparam TE Tensor expression type
 * @tparam D Dimension
 * @tparam R Tensor rank
 * @tparam dataType data type of tensor entries
 * @tparam assignable flag to indicate whether underlying expression is assignable
 */
template<typename TE, int D, int R, typename dataType, bool assignable>
class TensorExpression
{
public:
    //Sequence holding strides of each index through memqory
    using stride_sequence = typename geometric_sequence<R, D, 1>::type;

    // Compute required storage for a tensor of dimension D and rank R
    static constexpr int tensor_storage(int N)
    {
        return (N == 0) ? 1 : D * tensor_storage(N - 1);
    }


    int dim() const
    {
        return D;
    }

    int rank() const
    {
        return R;
    }


    // Down-casting functions
    TE const &self() const
    {
        return static_cast<TE const &>(*this);
    }

    TE &self()
    {
        return static_cast<TE &>(*this);
    }


    // Compute linear index of a (i,j,k...) component
    template<int S, typename INT>
    inline int index(sequence<S>, INT i) const noexcept
    {
        return i * S;
    }

    template<int S, int ...SS, typename INT, typename ...Args>
    inline int index(sequence<S, SS...>, INT i, Args ...args) const noexcept
    {
        return S * i + index(sequence<SS...>(), args...);
    }


    // Access tensor via linear indexing
    inline dataType linearIndexing(int idx) const noexcept
    {
        return self().linearIndexing(idx);
    }

    template<bool dummy_bool = assignable, typename = typename std::enable_if<dummy_bool>::type>
    inline dataType &linearIndexing(int idx) noexcept
    {
        return self().linearIndexing(idx);
    }


    // Parenthesis indexing operator
    template<typename ...Args,
            typename=typename std::enable_if<all_of_type<int>(type_list<Args...>())>::type>
    dataType operator()(Args... args) const noexcept { return linearIndexing(index(stride_sequence(),args...)); }


    template<bool dummy_bool = assignable, typename = typename std::enable_if<dummy_bool>::type,
            typename ...Args,
            typename=typename std::enable_if<all_of_type<int>(type_list<Args...>())>::type>
    dataType &operator()(Args... args) noexcept { return linearIndexing(index(stride_sequence(), args...)); }


    //Slicing indexing operator
    template<typename ...Args,
            typename=typename std::enable_if<contains<slice_sentinel>(type_list<Args...>())>::type
    >
    auto operator()(Args ...args) const noexcept
    {
        return make_slice(*this,
                          make_slice_offset(stride_sequence(), args...),
                          make_slice_strides(stride_sequence(), args...));
    }


    // Iterator classes for looping over all tensor elements in memory order
    class TensorExpressionIterator
    {
    public:
        TE &te_ref;
        int offset;

        TensorExpressionIterator(TensorExpression<TE, D, R, dataType, assignable> &te, int off)
                : te_ref(te.self()), offset(off) {}

        inline bool operator!=(const TensorExpressionIterator &other)
        {
            return offset != other.offset;
        }

        inline auto operator++()
        {
            ++offset;
            return *this;
        }


        inline dataType operator*() const noexcept
        {
            return te_ref.linearIndexing(offset);
        }

        template<bool dummy_bool = assignable, typename = typename std::enable_if<dummy_bool>::type>
        inline dataType &operator*() noexcept
        {
            return te_ref.linearIndexing(offset);
        }
    };


    class ConstTensorExpressionIterator
    {
    public:
        const TE &te_ref;
        int offset;

        ConstTensorExpressionIterator(const TensorExpression<TE, D, R, dataType, assignable> &te, int off)
                : te_ref(te.self()), offset(off) {}

        inline bool operator!=(const TensorExpressionIterator &other) const noexcept
        {
            return offset != other.offset;
        }

        inline auto operator++()
        {
            ++offset;
            return *this;
        }


        dataType operator*() const noexcept
        {
            return te_ref.linearIndexing(offset);
        }

        template<bool dummy_bool = assignable, typename = typename std::enable_if<dummy_bool>::type>
        dataType const &operator*() const noexcept
        {
            return te_ref.linearIndexing(offset);
        }
    };

    auto begin() const noexcept
    {
        return ConstTensorExpressionIterator(*this, 0);
    }

    auto end() const noexcept
    {
        return ConstTensorExpressionIterator(*this, tensor_storage(R));
    }

    auto begin() noexcept
    {
        return TensorExpressionIterator(*this, 0);
    }

    auto end() noexcept
    {
        return TensorExpressionIterator(*this, tensor_storage(R));
    }

};


template<template<typename, int, int, typename, bool> class TE, typename T, int D, int R,
        typename dataType, bool assignable, int ...PARENT_STRIDES>
auto make_slice(const TE<T, D, R, dataType, assignable> &te, int offset, sequence<PARENT_STRIDES...>)
{
    return TensorSlice<T, D, sizeof...(PARENT_STRIDES), dataType,
            false, PARENT_STRIDES...>(te, offset, sequence<PARENT_STRIDES...>());

};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSOREXPRESSION_HPP
