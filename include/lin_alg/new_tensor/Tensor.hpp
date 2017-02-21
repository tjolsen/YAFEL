//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_TENSOR_HPP
#define YAFEL_TENSOR_HPP

#include "TensorExpression.hpp"
#include "mp_utils/sequences.hpp"
#include <array>

YAFEL_NAMESPACE_OPEN

template<int D, int R, typename dataType=double>
class Tensor : public TensorExpression<Tensor<D, R, dataType>, D, R, dataType>
{

public:
    //Sequence holding strides of each index through memory
    using stride_sequence = typename geometric_sequence<R,D,1>::type;

    template<typename ...Args>
    dataType operator()(Args... args) const
    {
        static_assert(sizeof...(Args) == R, "Rank/argument mismatch");
        return data[index(stride_sequence(), args...)];
    }

    template<typename ...Args>
    dataType &operator()(Args... args)
    {
        static_assert(sizeof...(Args) == R, "Rank/argument mismatch");
        return data[index(stride_sequence(), args...)];
    }

    static constexpr int tensor_storage(int N) { return (N==0) ? 1 : D*tensor_storage(N-1); }

    template<int S, typename INT>
    int index(sequence<S>, INT i) { return i*S; }

    template<int S, int ...SS, typename INT, typename ...Args>
    int index(sequence<S, SS...>, INT i, Args ...args) { return S*i + index(sequence<SS...>(), args...); }


    // array to store tensor data
    std::array<dataType, tensor_storage(R)> data;
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSOR_HPP
