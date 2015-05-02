#ifndef _YAFEL_QUADRATURERULE_HPP
#define _YAFEL_QUADRATURERULE_HPP

#include "yafel_globals.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class QuadratureRule {

protected:
  std::vector<double> nodes;
  std::vector<double> weights;

public:
  QuadratureRule(unsigned polynomialOrder);

  inline unsigned get_nqp() const {return nodes.size();}
  inline double node(unsigned i) {return nodes[i];}
  inline double weight(unsigned i) {return weights[i];}
};


YAFEL_NAMESPACE_CLOSE

#endif
