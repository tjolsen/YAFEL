#ifndef _YAFEL_GAUSSLEGENDREQUADRATURE_HPP
#define _YAFEL_GAUSSLEGENDREQUADRATURE_HPP

#include "yafel_globals.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class GaussLegendreQuadrature {
  
private:
  std::vector<double> nodes;
  std::vector<double> weights;
  
public:
  GaussLegendreQuadrature(unsigned polynomialOrder);
  
  inline unsigned get_nqp() const { return nodes.size(); }
  inline double node(unsigned i) const {return nodes[i];};
  inline double weight(unsigned i) const {return weights[i];};
};

YAFEL_NAMESPACE_CLOSE

#endif
