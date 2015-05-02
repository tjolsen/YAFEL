#ifndef _YAFEL_GAUSSLOBATTOQUADRATURE_HPP
#define _YAFEL_GAUSSLOBATTOQUADRATURE_HPP

#include "yafel_globals.hpp"
#include "utils/QuadratureRule.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class GaussLobattoQuadrature : public QuadratureRule {
  
public:
  GaussLobattoQuadrature(unsigned Npoints);
  
};

YAFEL_NAMESPACE_CLOSE

#endif
