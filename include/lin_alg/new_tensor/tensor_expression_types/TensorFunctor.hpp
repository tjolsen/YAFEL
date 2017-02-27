//
// Created by tyler on 2/26/17.
//

#ifndef YAFEL_TENSORFUNCTOR_HPP
#define YAFEL_TENSORFUNCTOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/TensorExpression.hpp"
#include "lin_alg/new_tensor/mp_utils/sequences.hpp"
#include "lin_alg/new_tensor/mp_utils/mp_tuple_utils.hpp"
#include <utility>

YAFEL_NAMESPACE_OPEN

template<typename LAMBDA, int D, int R, typename dt=double>
class TensorFunctor : public TensorExpression<TensorFunctor<LAMBDA, D, R, dt>, D, R, dt, false>
{
public:
    using super = TensorExpression<TensorFunctor<LAMBDA, D, R, dt>, D, R, dt, false>;
    using idx_tuple = typename homogeneous_tuple<int, R>::type;

    LAMBDA function;

    TensorFunctor(LAMBDA &&func) : function(func) {}

    // easier to call via operator(), so will provide that,
    // then provide method to convert from linear index to argument pack
    template<typename ...Args>
    inline dt operator()(Args... args) const noexcept
    {
        static_assert(sizeof...(Args) == R, "Error: Incorrect number of arguments");
        return function(std::forward<Args...>(args)...);
    }

    inline dt linearIndexing(int idx) const noexcept
    {
        idx_tuple args;
        linear_to_index(args, sequence<0>(), typename super::stride_sequence(), idx);
        return apply_args(args, typename counting_sequence<R>::type());
    }

    template<int ...SS>
    inline dt apply_args(const idx_tuple &args, sequence<SS...>) const noexcept
    {
        return function(std::get<SS>(args)...);
    }

    template<int I, int S, int ...SS>
    inline void linear_to_index(idx_tuple &args, sequence<I>, sequence<S, SS...>, int idx) const noexcept
    {
        std::get<I>(args) = idx / S;
        linear_to_index(args, sequence<I + 1>(), sequence<SS...>(), idx % S);
    }

    template<int I>
    inline void linear_to_index(idx_tuple &, sequence<I>, sequence<>, int) const noexcept {}
};

template<typename LAMBDA, int D, int R, typename dt=double>
auto make_TensorFunctor(LAMBDA &&lambda, sequence<D, R>, type_list<dt>)
{
    return TensorFunctor<LAMBDA, D, R, dt>(std::forward<LAMBDA>(lambda));
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORFUNCTOR_HPP
