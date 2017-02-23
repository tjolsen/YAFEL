//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_TENSORSLICE_HPP
#define YAFEL_TENSORSLICE_HPP

#include "lin_alg/new_tensor/TensorExpression.hpp"

YAFEL_NAMESPACE_OPEN

/**
 * \class TensorSlice
 * \brief Slice a tensor akin to Matlab's C(i,j,:,:) notation.
 *
 * Create a lower-rank slice of a tensor, retaining the "assignability"
 * of the underlying TensorExpression. This slice should then be iterable
 * using standard begin(), end() functions (ie, a TensorSlice is a TensorExpression).
 */
template<typename TE, int D, int R, typename dataType, bool assignable, int ...PARENT_STRIDES>
class TensorSlice : public TensorExpression<TensorSlice<TE,D,R,dataType,assignable,PARENT_STRIDES...>, D,R,dataType,assignable>
{
public:
    using super = TensorExpression<TensorSlice<TE,D,R,dataType,assignable,PARENT_STRIDES...>, D,R,dataType,assignable>;
    using super::operator=;

    using parent_strides = sequence<PARENT_STRIDES...>;



    TE &te_ref;
    int offset;

    template<int TR, bool Tb>
    TensorSlice(TensorExpression<TE,D,TR,dataType,Tb> &te, int offset, parent_strides)
            : te_ref(te.self()), offset(offset)
    {
        static_assert(sizeof...(PARENT_STRIDES) == R, "PARENT_STRIDES size must equal rank");
    }

    template<typename TT, typename dt2, bool Tb>
    auto& operator=(const TensorExpression<TT,D,R,dt2,Tb> &rhs) noexcept
    {
        auto rit=rhs.begin();
        for(auto it=super::begin(); it!=super::end(); ++it, ++rit) {
            *it = *rit;
        }

        return *this;
    };


    inline dataType linearIndexing(int idx) const noexcept
    {
        return te_ref.linearIndexing(offset+local_to_parent(idx,typename super::stride_sequence(), parent_strides()));
    }

    template<bool dummy_bool=assignable, typename=typename std::enable_if<dummy_bool>::type>
    inline dataType &linearIndexing(int idx) noexcept
    {
        return te_ref.linearIndexing(offset+local_to_parent(idx,typename super::stride_sequence(), parent_strides()));
    }

    inline int local_to_parent(int, sequence<>, sequence<>) const noexcept { return 0; }

    template<int S, int ...SS, int P, int ...PP>
    inline int local_to_parent(int idx, sequence<S,SS...>, sequence<P, PP...>) const noexcept
    {
        return P*(idx/S) + local_to_parent(idx%S, sequence<SS...>(), sequence<PP...>());
    }
};




template<typename TE, int D, int R, typename dataType, bool assignable, int ...PARENT_STRIDES>
class ConstTensorSlice : public TensorExpression<ConstTensorSlice<TE,D,R,dataType,assignable,PARENT_STRIDES...>, D,R,dataType,assignable>
{
public:
    using super = TensorExpression<ConstTensorSlice<TE,D,R,dataType,assignable,PARENT_STRIDES...>, D,R,dataType,assignable>;
    using parent_strides = sequence<PARENT_STRIDES...>;

    const TE &te_ref;
    int offset;

    template<int TR, bool Tb>
    ConstTensorSlice(const TensorExpression<TE,D,TR,dataType,Tb> &te, int offset, parent_strides)
            : te_ref(te.self()), offset(offset)
    {
        static_assert(sizeof...(PARENT_STRIDES) == R, "PARENT_STRIDES size must equal rank");
    }


    inline dataType linearIndexing(int idx) const noexcept
    {
        return te_ref.linearIndexing(offset+local_to_parent(idx,typename super::stride_sequence(), parent_strides()));
    }

    /*
    template<bool dummy_bool=assignable, typename=typename std::enable_if<dummy_bool>::type>
    inline dataType &linearIndexing(int idx) noexcept
    {
        return te_ref.linearIndexing(offset+local_to_parent(idx,typename super::stride_sequence(), parent_strides()));
    }
    */

    inline int local_to_parent(int, sequence<>, sequence<>) const noexcept { return 0; }

    template<int S, int ...SS, int P, int ...PP>
    inline int local_to_parent(int idx, sequence<S,SS...>, sequence<P, PP...>) const noexcept
    {
        return P*(idx/S) + local_to_parent(idx%S, sequence<SS...>(), sequence<PP...>());
    }
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORSLICE_HPP
