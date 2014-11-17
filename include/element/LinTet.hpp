#ifndef _YAFEL_LINTET_HPP
#define _YAFEL_LINTET_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

class LinTet : public Element {

public:
  LinTet(int dofpn);
  
  double shape_value_xi(int node, const Vector &xi) const;
  double shape_grad_xi(int node, int component, const Vector &xi) const;

};

YAFEL_NAMESPACE_CLOSE

#endif
