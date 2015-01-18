#ifndef _YAFEL_LINTET_HPP
#define _YAFEL_LINTET_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

class LinTet : public Element {

public:
  LinTet(const DoFManager &dofm);
  
  double shape_value_xi(unsigned node, const Vector &xi) const;
  double shape_grad_xi(unsigned node, unsigned component, const Vector &xi) const;

};

YAFEL_NAMESPACE_CLOSE

#endif
