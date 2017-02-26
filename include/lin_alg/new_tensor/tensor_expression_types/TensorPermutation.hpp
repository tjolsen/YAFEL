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
        static_assert(sizeof...(IDX_PERM) == R, "Invalid permutation length");
    }

    template<int S, int ...SS, int I, int ...II>
    int local_to_parent(int idx, stride_sequence, sequence<S, SS...>, sequence<I, II...>) const noexcept
    {
        return index_at(stride_sequence(), I)*(idx/S) + local_to_parent(idx%S, stride_sequence(),
                                                                        sequence<SS...>(), sequence<II...>());
    };

    int local_to_parent(int, stride_sequence, sequence<>, sequence<>) const noexcept
    {
        return 0;
    }

    inline dt linearIndexing(int idx) const noexcept
    {
        return te_ref.linearIndexing(local_to_parent(idx,stride_sequence(), stride_sequence(), idx_perm()));
    }


};


template<typename TE, int D, int R, typename dt, bool b, int ...IDX_PERM>
auto permute(const TensorExpression<TE,D,R,dt,b>& te, sequence<IDX_PERM...>) {
    return TensorPermutation<TE,D,R,dt,b,IDX_PERM...>(te,sequence<IDX_PERM...>());
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORPERMUTATION_HPP
