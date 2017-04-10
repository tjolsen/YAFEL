//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_TENSORSCALED_HPP
#define YAFEL_TENSORSCALED_HPP

#include "lin_alg/tensor/TensorExpression.hpp"

YAFEL_NAMESPACE_OPEN

template<typename TE, typename scaleT, int D, int R, typename dataType>
class TensorScaled : public TensorExpression<TensorScaled<TE,scaleT,D,R,dataType>,D,R,decltype(dataType()*scaleT()),false>
{
public:
    using super = TensorExpression<TensorScaled<TE,scaleT,D,R,dataType>,D,R,decltype(dataType()*scaleT()),false>;
    const TE& te_ref;
    const scaleT scale;

    template<typename dt1, typename dt2, bool b>
    TensorScaled(const TensorExpression<TE,D,R,dt1,b> &te, dt2 s)
            : te_ref(te.self()), scale(s)
    {}


    inline dataType linearIndexing(int idx) const noexcept
    {
        return scale*te_ref.linearIndexing(idx);
    }
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORSCALED_HPP
