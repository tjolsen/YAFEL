#ifndef _YAFEL_ELEMENTFACTORY_HPP
#define _YAFEL_ELEMENTFACTORY_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "element/LinHex.hpp"
#include "element/LinQuad.hpp"
#include "element/LinTri.hpp"
#include "element/NullElement.hpp"
#include "element/LinTet.hpp"
#include "element/LinLine.hpp"
//#include "element/LagrangeLine.hpp"
#include "mesh/GenericMesh.hpp"
#include "utils/DoFManager.hpp"
#include "utils/ElementType.hpp"

YAFEL_NAMESPACE_OPEN

template<typename MTYPE, unsigned NSD>
class ElementFactory {

public:
    ElementFactory()=delete;
    ElementFactory(const GenericMesh<MTYPE,NSD> &M, const DoFManager &dofm);

    Element<NSD> & getElement(size_type elnum); //return NullElement if not implemented

    inline const GenericMesh<MTYPE,NSD> & getMesh() const {return M;}
    inline size_type n_dof() const { return dof_per_node()*(M.n_nodes()); }
    inline size_type dof_per_node() const { return _dof_per_node; }


private:
    const GenericMesh<MTYPE,NSD> &M;
    LinLine<NSD> linear_line;
    LinQuad<NSD> linear_quad;
    LinTri<NSD> linear_tri;
    LinHex<NSD> linear_hex;
    LinTet<NSD> linear_tet;
    NullElement<NSD> null_element;
    const DoFManager & DOFM;
    size_type _dof_per_node;
};


/*
 * Implementation starts here
 */

template<typename MTYPE, unsigned NSD>
ElementFactory<MTYPE,NSD>::ElementFactory(const GenericMesh<MTYPE,NSD> &m, const DoFManager &dofm)
    : M(m),
      linear_line(dofm),
      linear_quad(dofm),
      linear_tri(dofm),
      linear_hex(dofm),
      linear_tet(dofm),
      null_element(dofm),
      DOFM(dofm),
      _dof_per_node(dofm.dof_per_node())
{

    linear_line.init_element();
    linear_quad.init_element();
    linear_tri.init_element();
    linear_hex.init_element();
    linear_tet.init_element();
}

template<typename MTYPE, unsigned NSD>
Element<NSD> & ElementFactory<MTYPE,NSD>::getElement(size_type elnum) {
  
    ElementType t = M.element_type(elnum);
  
    switch(t) {
    case ElementType::LINEAR_LINE:
	return linear_line;
    case ElementType::LINEAR_QUAD:
        return linear_quad;
    case ElementType::LINEAR_TRI:
        return linear_tri;
    case ElementType::LINEAR_HEX:
        return linear_hex;
    case ElementType::LINEAR_TET:
        return linear_tet;
    default:
        return null_element;
    }
}

YAFEL_NAMESPACE_CLOSE

#endif //#ifndef _YAFEL_ELEMENTFACTORY_HPP
