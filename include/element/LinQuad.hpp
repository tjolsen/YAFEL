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
template<unsigned NSD, typename std::enable_if<NSD>=2>::type >
LinQuad<NSD>::LinQuad(const DoFManager &dofm) : Element(dofm, 2, 4, 4*dofm.getDofPerNode(), 9, 4)
{

  //set up parent element xi_0 vectors
  xi_0.clear();
  xi_0.emplace_back(coordinate_type{-1, -1});
  xi_0.emplace_back(coordinate_type{ 1, -1});
  xi_0.emplace_back(coordinate_type{ 1,  1});
  xi_0.emplace_back(coordinate_type{-1,  1});

  //assign gauss points and weights
  quad_points.clear();
  gauss_weights.clear();
  double a = 1.0/sqrt(3.0);
  v(0) = -a; v(1) = -a; quad_points.push_back(v); gauss_weights.push_back(1.0);
  v(0) = a; v(1) = -a; quad_points.push_back(v); gauss_weights.push_back(1.0);
  v(0) = a; v(1) = a; quad_points.push_back(v); gauss_weights.push_back(1.0);
  v(0) = -a; v(1) = a; quad_points.push_back(v); gauss_weights.push_back(1.0);
}


template<unsigned NSD>
double LinQuad<NSD>::shape_value_xi(size_type node, const coordinate_type &xi) const {
  
  double val = 1.0;
  for(unsigned i=0; i<NSD; ++i) {
    val *= (1.0/2.0) * (xi_0[node](i)*xi(i) + 1);
  }
  
  return val;
}

template<unsigned NSD>
double LinQuad::shape_grad_xi(size_type node, size_type comp, const Vector &xi) const {
  
  double val = 1.0;
  for(size_type i=0; i<NSD; ++i) {
    if(comp == i)
      continue;
    val *= (1.0/2.0) * (xi_0[node](i)*xi(i)+1);
  }
  
  return val * xi_0[node](comp)/(2.0);
}


YAFEL_NAMESPACE_CLOSE

#endif
