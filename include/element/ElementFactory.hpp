#ifndef __YAFEL_ELEMENTFACTORY_HPP
#define __YAFEL_ELEMENTFACTORY_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "element/LinQuad.hpp"
#include "element/LinTri.hpp"
#include "element/LinTet.hpp"
#include "element/LinLine.hpp"
#include "element/LagrangeLine.hpp"
#include "mesh/GenericMesh.hpp"
#include "utils/DoFManager.hpp"
#include "utils/ElementType.hpp"

YAFEL_NAMESPACE_OPEN

template<typename MTYPE, unsigned NSD>
class ElementFactory {

public:
  using size_type = typename GenericMesh<MTYPE,NSD>::size_type;

  ElementFactory()=delete;
  ElementFactory(const GenericMesh<MTYPE,NSD> &M, const DoFManager &dofm);

  Element & getElement(size_type elnum); //return NullElement if not implemented

  inline const GenericMesh<MTYPE,NSD> & getMesh() const {return M;}
  inline size_type n_dof() const { return dof_per_node*(M.n_nodes()); }
  inline size_type dof_per_node() const { return _dof_per_node; }


private:
  const GenericMesh<MTYPE,NSD> &M;
  LinQuad<NSD> linear_quad;
  LinTri<NSD> linear_tri;
  NullElement<NSD> null_element;
  DoFManager DOFM;
  size_type _dof_per_node;
};


/*
 * Implementation starts here
 */

template<typename MTYPE, unsigned NSD>
ElementFactory(const GenericMesh<MTYPE,NSD> &M, const DoFManager &dofm)
  : linear_quad(dofm),
    linear_tri(dofm),
    null_element(dofm),
    DOFM(dofm),
    dof_per_node(dofm.dof_per_node())
{}

template<typename MTYPE, unsigned NSD>
Element & getElement(size_type elnum) {
  
  ElementType t = M.element_type(elnum);

  switch(t) {
  case ElementType::LINEAR_QUAD:
    return linear_quad;
  default:
    return null_element;
  }
}

YAFEL_NAMESPACE_CLOSE

#endif //#ifndef _YAFEL_ELEMENTFACTORY_HPP
