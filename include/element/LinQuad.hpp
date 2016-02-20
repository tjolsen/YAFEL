#ifndef _YAFEL_LINQUAD_HPP
#define _YAFEL_LINQUAD_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"
#include <vector>
#include <type_traits>

YAFEL_NAMESPACE_OPEN

template<NSD>
class LinQuad : public Element<NSD> {

public:
  using size_type = typename Element<NSD>::size_type;
  using coordinate_type = typename Element<NSD>::coordinate_type;


  LinQuad(const DoFManager &dofm);

  double shape_value_xi(size_type node, const coordinate_type &xi) const;
  double shape_grad_xi(size_type node, size_type component, const coordinate_type &xi) const;
  
};


/*
 * LinQuad Element Implementation
 */
template<unsigned NSD, std::enable_if<NSD>=2> >
LinQuad<NSD>::LinQuad(const DoFManager &dofm) : Element(dofm, 2, 4, 4*dofm.getDofPerNode(), 9, 4)
{

  coordinate_type v;
  
  //set up parent element xi_0 vectors
  xi_0.clear();
  v(0) = -1; v(1) = -1; xi_0.push_back(v);
  v(0) = 1; v(1) = -1; xi_0.push_back(v);
  v(0) = 1; v(1) = 1; xi_0.push_back(v);
  v(0) = -1; v(1) = 1; xi_0.push_back(v);

  //assign gauss points and weights
  quad_points.clear();
  gauss_weights.clear();
  double a = 1.0/sqrt(3.0);
  v(0) = -a; v(1) = -a; quad_points.push_back(v); gauss_weights.push_back(1.0);
  v(0) = a; v(1) = -a; quad_points.push_back(v); gauss_weights.push_back(1.0);
  v(0) = a; v(1) = a; quad_points.push_back(v); gauss_weights.push_back(1.0);
  v(0) = -a; v(1) = a; quad_points.push_back(v); gauss_weights.push_back(1.0);
}


YAFEL_NAMESPACE_CLOSE

#endif
