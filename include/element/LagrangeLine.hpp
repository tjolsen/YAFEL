#ifndef _YAFEL_LAGRANGELINE_HPP
#define _YAFEL_LAGRANGELINE_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

class LagrangeLine : public Element {


public:
  LagrangeLine(const DoFManager &dofm, unsigned polynomialOrder);
  
  double shape_value_xi(unsigned node, const Vector &xi) const;
  double shape_grad_xi(unsigned node, unsigned component, const Vector &xi) const;

};


YAFEL_NAMESPACE_CLOSE

#endif
