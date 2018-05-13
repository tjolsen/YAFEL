//
// Created by tyler on 3/9/17.
//

#ifndef YAFEL_UNARYFUNCTIONS_HPP
#define YAFEL_UNARYFUNCTIONS_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorCwiseUnaryOp.hpp"
#include "lin_alg/tensor/tensor_expression_types/UnaryOperations.hpp"

YAFEL_NAMESPACE_OPEN

template<typename T, int D, int R, typename dt, bool b>
auto sqrt(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,Sqrt>(te);
}


template<typename T, int D, int R, typename dt, bool b>
auto sign(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,Sign>(te);
}


template<typename T, int D, int R, typename dt, bool b>
auto sin(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,Sin>(te);
}


template<typename T, int D, int R, typename dt, bool b>
auto cos(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,Cos>(te);
}

template<typename T, int D, int R, typename dt, bool b>
auto round(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,Round>(te);
}

template<typename T, int D, int R, typename dt, bool b>
auto floor(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,Floor>(te);
}

template<typename T, int D, int R, typename dt, bool b>
auto ceil(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,Ceil>(te);
}

template<typename T, int D, int R, typename dt, bool b>
auto abs(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,Abs>(te);
};

template<typename T, int D, int R, typename dt, bool b>
auto isnan(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,IsNan>(te);
};

template<typename CastType, typename T, int D, int R, typename dt, bool b>
auto tensor_cast(const TensorExpression<T,D,R,dt,b> &te) noexcept
{
    return TensorCwiseUnaryOp<T,D,R,dt,b,typename TensorCaster<CastType>::Cast>(te);
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_UNARYFUNCTIONS_HPP
