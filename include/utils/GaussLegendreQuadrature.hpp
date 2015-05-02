#ifndef _YAFEL_GAUSSLEGENDREQUADRATURE_HPP
#define _YAFEL_GAUSSLEGENDREQUADRATURE_HPP

#include "yafel_globals.hpp"
#include "utils/QuadratureRule.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class GaussLegendreQuadrature : public QuadratureRule {
  
public:
  GaussLegendreQuadrature(unsigned polynomialOrder);
  
};

YAFEL_NAMESPACE_CLOSE

#endif
