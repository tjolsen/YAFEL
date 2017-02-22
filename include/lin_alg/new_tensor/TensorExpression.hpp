//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_TENSOREXPRESSION_HPP
#define YAFEL_TENSOREXPRESSION_HPP

#include "yafel_globals.hpp"
#include "mp_utils/sequences.hpp"
#include "mp_utils/sequence_functions.hpp"

#include <type_traits>

YAFEL_NAMESPACE_OPEN

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
    int index(sequence<S>, INT i)
    {
        return i * S;
    }

    template<int S, int ...SS, typename INT, typename ...Args>
    int index(sequence<S, SS...>, INT i, Args ...args)
    {
        return S * i + index(sequence<SS...>(), args...);
    }


    // Access tensor via linear indexing
    inline dataType linearIndexing(int idx) const
    {
        return self().linearIndexing(idx);
    }

    template<bool dummy_bool = assignable, typename = typename std::enable_if<dummy_bool>::type>
    inline dataType &linearIndexing(int idx)
    {
        return self().linearIndexing(idx);
    }


    // Parenthesis indexing operator
    template<typename ...Args>
    dataType operator()(Args... args) const
    { return linearIndexing(stride_sequence(), index(args...)); }

    template<bool dummy_bool = assignable, typename = typename std::enable_if<dummy_bool>::type,
            typename ...Args>
    dataType &operator()(Args... args)
    { return linearIndexing(index(stride_sequence(), args...)); }


    // Iterator classes for looping over all tensor elements in memory order
    class TensorExpressionIterator
    {
    public:
        TE &te_ref;
        int offset;

        TensorExpressionIterator(TensorExpression<TE, D, R, dataType, assignable> &te, int off)
                : te_ref(te.self()), offset(off)
        {}

        bool operator!=(const TensorExpressionIterator &other)
        {
            return offset != other.offset;
        }

        auto operator++()
        {
            ++offset;
            return *this;
        }


        dataType operator*() const
        {
            return te_ref.linearIndexing(offset);
        }

        template<bool dummy_bool = assignable, typename = typename std::enable_if<dummy_bool>::type>
        dataType &operator*()
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
                : te_ref(te.self()), offset(off)
        {}

        bool operator!=(const TensorExpressionIterator &other)
        {
            return offset != other.offset;
        }

        auto operator++()
        {
            ++offset;
            return *this;
        }


        dataType operator*() const
        {
            return te_ref.linearIndexing(offset);
        }

        template<bool dummy_bool = assignable, typename = typename std::enable_if<dummy_bool>::type>
        dataType const &operator*()
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


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSOREXPRESSION_HPP
