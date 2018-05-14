//
// Created by tyler on 5/14/18.
//

#ifndef YAFEL_TENSORCONSTANT_HPP
#define YAFEL_TENSORCONSTANT_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"

YAFEL_NAMESPACE_OPEN

/**
 * \class TensorConstant
 * \brief A tensor of dimension D and rank R filled with a single value.
 *
 * This is intended to represent a constant tensor without incurring the storage
 * associated with a full Tensor.
 * @tparam D tensor dimension
 * @tparam R tensor rank
 * @tparam dt datatype
 */
template<int D, int R, typename dt = double>
class TensorConstant : public TensorExpression<TensorConstant<D,R,dt>, D, R, dt, false>
{
public:
    template<typename T=dt, typename=typename std::enable_if_t<std::is_convertible_v<T,dt>>>
    TensorConstant(T const& t = dt{0}) : value(t) {}

    YAFEL_ALWAYS_INLINE dt linearIndexing(int) const noexcept { return value; }
private:
    dt value;
};

/**
 * Useful variable template to represent a "zero" tensor, without the storage of an actual Tensor
 */
template<int D, int R, typename dt = double>
TensorConstant<D,R,dt> TensorZero = TensorConstant<D,R,dt>(dt{0});

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORCONSTANT_HPP
