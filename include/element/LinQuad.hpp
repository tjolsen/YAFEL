#ifndef _YAFEL_LINQUAD_HPP
#define _YAFEL_LINQUAD_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class LinQuad : public Element {

public:
  LinQuad(unsigned dofpn);

  double shape_value_xi(unsigned node, const Vector &xi) const;
  double shape_grad_xi(unsigned node, unsigned component, const Vector &xi) const;
  
};

YAFEL_NAMESPACE_CLOSE

#endif
