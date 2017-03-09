//
// Created by tyler on 2/21/17.
//

#include "lin_alg/new_tensor/Tensor.hpp"
#include "lin_alg/new_tensor/tensor_functions/operators.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSlice.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorContraction.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorFunctor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorMap.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorDyadicProduct.hpp"
#include "lin_alg/new_tensor/tensor_functions/tensor_dot.hpp"
#include "lin_alg/new_tensor/tensor_functions/UnaryFunctions.hpp"
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"
#include "lin_alg/new_tensor/mp_utils/slice_mp_utils.hpp"
#include "utils/DualNumber.hpp"


#include <iostream>
#include <cmath>
#include <vector>

using namespace yafel;
//using namespace std;


struct NewNum
{
    int a, b;

    NewNum() : NewNum(0,0) {}

    template<typename T>
    NewNum(T t, T u=0) : a(t), b(u) {}

    NewNum operator*(const NewNum &rhs)
    {
        return NewNum{a*rhs.a, b*rhs.b};
    }

    NewNum operator-() const {return NewNum(-a, -b); }

    NewNum &operator=(int rhs) {a=rhs; b=rhs; return *this;}

    NewNum &operator+=(const NewNum &rhs) { a+= rhs.a; b += rhs.b; return *this; }

    NewNum operator+(const NewNum &rhs) {NewNum tmp(*this); return tmp += rhs;}
};


template<>
struct ScalarTraits<NewNum>
{
    static constexpr bool isYafelScalar() { return true; }
};



int main()
{

    //Tensor<4,1,DualNumber<double>> x,y;
    Tensor<4,1,NewNum> x,y;
    //Tensor<4,1,double> x,y;
    int count = 1;
    auto yit = y.begin();
    for(auto &xi : x) {
        *yit = xi = count++;
        ++yit;
    }

    //auto z = sqrt(x).eval();

    //for(auto zi : z) {
    //std::cout << zi << std::endl;
    //}

    std::cout << ScalarTraits<NewNum>::isYafelScalar() << std::endl;

    auto yy = dot(y,-y);

    return yy.a>0;
}
