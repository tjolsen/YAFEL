//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_TENSOR_HPP
#define YAFEL_TENSOR_HPP

#include "lin_alg/tensor/TensorExpression.hpp"
#include <array>

YAFEL_NAMESPACE_OPEN

template<int D, int R, typename dataType=double>
class Tensor : public TensorExpression<Tensor<D, R, dataType>, D, R, dataType, true>
{

public:
    // super-type
    using super = TensorExpression<Tensor<D, R, dataType>, D, R, dataType, true>;

    // Initialize a tensor with all elements set to "d"
    template<typename dt,
            typename=typename std::enable_if<std::is_convertible<dataType, dt>::value>::type>
    Tensor(dt d = dt(0))
    {
        for (int i = 0; i < super::tensor_storage(R); ++i) {
            data[i] = d;
        }
    }

    // Initialize an empty tensor with all zeros
    Tensor() : Tensor(dataType(0)) {}

    template<typename ...DT,
            typename=typename std::enable_if<2 <= sizeof...(DT)>::type>
    Tensor(DT ...args) : data{{std::forward<dataType>(args)...}} {}


    // Construct from a TensorExpression
    template<typename TE, typename dt2, bool b, typename=typename std::enable_if<std::is_convertible<dataType, dt2>::value>::type>
    Tensor(const TensorExpression <TE, D, R, dt2, b> &rhs)
    {
        auto rit = rhs.begin();
        for (int i = 0; i < super::tensor_storage(R); ++i, ++rit) {
            data[i] = dataType(*rit);
        }
    }

    //implement indexing operations
    dataType linearIndexing(int idx) const noexcept
    {
        return data[idx];
    }

    dataType &linearIndexing(int idx) noexcept
    {
        return data[idx];
    }

    // array to store tensor data
    std::array<dataType, super::tensor_storage(R)> data;
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSOR_HPP
