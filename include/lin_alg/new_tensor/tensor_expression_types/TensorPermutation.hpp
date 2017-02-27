//
// Created by tyler on 2/26/17.
//

#ifndef YAFEL_TENSORPERMUTATION_HPP
#define YAFEL_TENSORPERMUTATION_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/TensorExpression.hpp"
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"

YAFEL_NAMESPACE_OPEN

/**
 * \class TensorPermutation
 * \brief Permute the indices of a tensor
 *
 * Permute the indices of a tensor expression, resulting in (eg) the
 * following behavior: result(i,j,k) := source(i,k,j)
 *
 * Any valid permutation of <0,1,...R> is permitted, and the resulting
 * TensorExpression is treated (for now) as NOT being assignable.
 *
 * @tparam TE TensorExpression type
 * @tparam D Tensor Dimension
 * @tparam R Tensor Rank (/order)
 * @tparam dt dataType
 * @tparam b bool flag indicating whether underlying TensorExpression is assignable
 */
template<typename TE, int D, int R, typename dt, bool b, int ...IDX_PERM>
class TensorPermutation : public TensorExpression<TensorPermutation<TE,D,R,dt,b,IDX_PERM...>,D,R,dt,false>
{
public:
    using super = TensorExpression<TensorPermutation<TE,D,R,dt,b,IDX_PERM...>,D,R,dt,false>;
    using stride_sequence = typename super::stride_sequence;
    using idx_perm = sequence<IDX_PERM...>;

    const TE& te_ref;

    TensorPermutation(const TensorExpression<TE,D,R,dt,b> &T, idx_perm)
            : te_ref(T.self())
    {
        validate_idx_perm();
    }

    template<int S, int ...SS, int I, int ...II>
    int local_to_parent(int idx, stride_sequence, sequence<S, SS...>, sequence<I, II...>) const noexcept
    {
        auto p_stride = index_at(stride_sequence(), sequence<index_at(idx_perm(), sequence<I>())>());
        return p_stride*(idx/S)
               + local_to_parent(idx%S, stride_sequence(), sequence<SS...>(), sequence<II...>());
    };

    int local_to_parent(int, stride_sequence, sequence<>, sequence<>) const noexcept
    {
        return 0;
    }

    inline dt linearIndexing(int idx) const noexcept
    {
        return te_ref.linearIndexing(local_to_parent(idx,stride_sequence(), stride_sequence(), typename counting_sequence<R>::type()));
    }


    void validate_idx_perm()
    {
        static_assert(sizeof...(IDX_PERM) == R, "Invalid permutation length");
        static_assert(all_ge(idx_perm(), sequence<0>()), "Error: IDX_PERM contains negative number");
        static_assert(all_lt(idx_perm(), sequence<R>()), "Error: IDX_PERM contains indices >= Rank. Remember, zero-based indexing!");
        static_assert(!contains_duplicates(idx_perm()), "Error: IDX_PERM contains duplicates. IDX_PERM must valid permutation of {0...R-1}.");
    }
};


template<typename TE, int D, int R, typename dt, bool b, int ...IDX_PERM>
auto permute(const TensorExpression<TE,D,R,dt,b>& te, sequence<IDX_PERM...>) {
    return TensorPermutation<TE,D,R,dt,b,IDX_PERM...>(te,sequence<IDX_PERM...>());
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORPERMUTATION_HPP
