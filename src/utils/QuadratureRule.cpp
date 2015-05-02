#include "utils/QuadratureRule.hpp"

YAFEL_NAMESPACE_OPEN

QuadratureRule::QuadratureRule(unsigned polynomialOrder) :
  nodes(polynomialOrder, 0.0), weights(polynomialOrder, 0.0)
{}

YAFEL_NAMESPACE_CLOSE
