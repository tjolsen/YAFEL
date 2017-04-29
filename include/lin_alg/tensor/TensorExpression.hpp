//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_TENSOREXPRESSION_HPP
#define YAFEL_TENSOREXPRESSION_HPP

#include "yafel_globals.hpp"
#include "utils/ScalarTraits.hpp"
#include "mp_utils/sequences.hpp"
#include "mp_utils/sequence_functions.hpp"
#include "mp_utils/TypeList.hpp"
#include "mp_utils/slice_mp_utils.hpp"

#include <type_traits>
#include <iterator>

YAFEL_NAMESPACE_OPEN

// Forward declaration of Tensor, TensorSlice and TensorPermutation so that they can be used from the base class
template<int D, int R, typename dataType>
class Tensor;

template<typename TE, int D, int R, typename dataType, bool assignable, int ...PARENT_STRIDES>
class TensorSlice;

template<typename TE, int D, int R, typename dataType, bool assignable, int ...PARENT_STRIDES>
class ConstTensorSlice;

template<typename TE, int D, int R, typename dt, bool b, int ...IDX_PERM>
class TensorPermutation;


template<template<typename, int, int, typename, bool> class TE, typename T, int D, int R,
        typename dataType, bool assignable, int ...PARENT_STRIDES, typename ...Args>
auto make_slice(const TE<T, D, R, dataType, assignable> &te, Args... args);

template<template<typename, int, int, typename, bool> class TE, typename T, int D, int R,
        typename dataType, bool assignable, int ...PARENT_STRIDES, typename ...Args>
auto make_const_slice(const TE<T, D, R, dataType, assignable> &te, Args... args);

template<template<typename, int, int, typename, bool> class TE, typename T, int D, int R,
        typename dt, int ...IDX_PERM>
auto permute(TE<T, D, R, dt, true> &te, sequence<IDX_PERM...>);

template<template<typename, int, int, typename, bool> class TE, typename T, int D, int R,
        typename dt, bool b, int ...IDX_PERM>
auto const_permute(const TE<T, D, R, dt, b> &te, sequence<IDX_PERM...>);


/**
 * \class TensorExpression
 * \brief Base class for tensor expressions. All others derive from this.
 * @tparam TE Tensor expression type
 * @tparam D Dimension
 * @tparam R Tensor rank
 * @tparam dataType data type of tensor entries
 * @tparam assignable flag to indicate whether underlying expression is assignable
 */
template<typename TE, int D, int R, typename dataType, bool assignable, typename=typename std::enable_if<ScalarTraits<dataType>::isYafelScalar()>::type>
class TensorExpression
{
public:
    //Sequence holding strides of each index through memqory
    using stride_sequence = typename geometric_sequence<R, D, 1>::type;

    // Compute required storage for a tensor of dimension D and rank R
    static constexpr int tensor_storage(int N)
    {
        int sz{1};
        while (N > 0) {
            sz *= D;
            --N;
        }
        return sz;
    }


    int dim() const { return D; }

    int rank() const { return R; }

    // Down-casting functions
    TE const &self() const { return static_cast<TE const &>(*this); }

    TE &self() { return static_cast<TE &>(*this); }

    // Evaluate the current expression into a Tensor<D,R,dataType>
    auto eval() const noexcept
    {
        return Tensor<D, R, dataType>(self());
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
    dataType operator()(Args... args) const noexcept { return linearIndexing(index(stride_sequence(), args...)); }


    template<bool dummy_bool = assignable, typename = typename std::enable_if<dummy_bool>::type,
            typename ...Args,
            typename=typename std::enable_if<all_of_type<int>(type_list<Args...>())>::type>
    dataType &operator()(Args... args) noexcept { return linearIndexing(index(stride_sequence(), args...)); }


    //Slicing indexing operator
    template<typename ...Args,
            typename=typename std::enable_if<contains<slice_sentinel>(type_list<Args...>())>::type
    >
    auto operator()(Args ...args) noexcept
    {
        return make_slice(*this,
                          make_slice_offset(stride_sequence(), args...),
                          make_slice_strides(stride_sequence(), args...));
    }

    template<typename ...Args,
            typename=typename std::enable_if<contains<slice_sentinel>(type_list<Args...>())>::type
    >
    auto operator()(Args ...args) const noexcept
    {
        return make_const_slice(*this,
                                make_slice_offset(stride_sequence(), args...),
                                make_slice_strides(stride_sequence(), args...));
    }

    /**
     * Permute the tensor expression
     * @return
     */
    template<int ...IDXS>
    auto perm() const noexcept
    {
        return const_permute(*this, sequence<IDXS...>());
    }

    template<int ...IDXS>
    auto perm() noexcept
    {
        return permute(*this, sequence<IDXS...>());
    }

    // Iterator classes for looping over all tensor elements in memory order
    class TensorExpressionIterator
    {
    public:
        using value_type = dataType;
        using reference = dataType&;

        TE &te_ref;
        int offset;

        TensorExpressionIterator(TensorExpression<TE, D, R, dataType, assignable> &te, int off)
                : te_ref(te.self()), offset(off) {}

        inline bool operator!=(const TensorExpressionIterator &other)
        {
            return offset != other.offset;
        }

        inline bool operator==(const TensorExpressionIterator &other)
        {
            return offset == other.offset;
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
        using value_type = dataType;
        using reference = const dataType&;
        using iterator_category = std::forward_iterator_tag;
        using pointer = const dataType*;
        using difference_type = int;

        const TE &te_ref;
        int offset;

        ConstTensorExpressionIterator(const TensorExpression<TE, D, R, dataType, assignable> &te, int off)
                : te_ref(te.self()), offset(off) {}

        inline bool operator!=(const ConstTensorExpressionIterator &other) const noexcept
        {
            return offset != other.offset;
        }
        inline bool operator==(const ConstTensorExpressionIterator &other) const noexcept
        {
            return offset == other.offset;
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

    ConstTensorExpressionIterator begin() const noexcept
    {
        return ConstTensorExpressionIterator(*this, 0);
    }

    ConstTensorExpressionIterator end() const noexcept
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


template<template<typename, int, int, typename, bool, typename> class TE, typename T, int D, int R,
        typename dataType, bool assignable, typename enabled, int ...PARENT_STRIDES>
auto make_slice(TE<T, D, R, dataType, assignable, enabled> &te, int offset, sequence<PARENT_STRIDES...>)
{
    return TensorSlice<T, D, sizeof...(PARENT_STRIDES), dataType,
            assignable, PARENT_STRIDES...>(te, offset, sequence<PARENT_STRIDES...>());

}


template<template<typename, int, int, typename, bool, typename> class TE, typename T, int D, int R,
        typename dataType, bool assignable, typename enabled, int ...PARENT_STRIDES>
auto make_const_slice(const TE<T, D, R, dataType, assignable, enabled> &te, int offset, sequence<PARENT_STRIDES...>)
{
    return ConstTensorSlice<T, D, sizeof...(PARENT_STRIDES), dataType,
            false, PARENT_STRIDES...>(te, offset, sequence<PARENT_STRIDES...>());

}

template<template<typename, int, int, typename, bool, typename> class TE, typename T, int D,
        int R, typename dt, typename enabled, int ...IDX_PERM>
auto permute(TE<T, D, R, dt, true, enabled> &te, sequence<IDX_PERM...>)
{
    return TensorPermutation<T, D, R, dt, true, IDX_PERM...>(te, sequence<IDX_PERM...>());
}

template<template<typename, int, int, typename, bool, typename> class TE, typename T, int D,
        int R, typename dt, bool b, typename enabled, int ...IDX_PERM>
auto const_permute(const TE<T, D, R, dt, b, enabled> &te, sequence<IDX_PERM...>)
{
    return TensorPermutation<T, D, R, dt, b, IDX_PERM...>(te, sequence<IDX_PERM...>());
}

YAFEL_NAMESPACE_CLOSE

/*
template<typename TETYPE>
struct std::iterator_traits<typename TETYPE::ConstTensorExpressionIterator>
{
    using iterType = typename TETYPE::ConstTensorExpressionIterator;
    using value_type = typename iterType::value_type;
    using reference = typename iterType::reference;
    using difference_type = int;

    using iterator_category = std::forward_iterator_tag;
};
*/
#endif //YAFEL_TENSOREXPRESSION_HPP
