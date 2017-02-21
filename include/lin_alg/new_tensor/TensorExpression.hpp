//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_TENSOREXPRESSION_HPP
#define YAFEL_TENSOREXPRESSION_HPP

#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

template<typename TE, int D, int R, typename dataType>
class TensorExpression
{
public:

    int dim() const
    { return D; }

    int rank() const
    { return R; }

    TE const &self() const
    { return static_cast<TE const &>(*this); }

    TE &self()
    { return static_cast<TE &>(*this); }

    template<typename ...Args>
    dataType operator()(Args... args) const
    { return self()(args...); }
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSOREXPRESSION_HPP
