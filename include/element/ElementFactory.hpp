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

YAFEL_NAMESPACE_OPEN

template<typename MTYPE, unsigned NSD>
class ElementFactory {

private:
  const GenericMesh<MTYPE,NSD> &M;
  LinQuad<NSD> *lq;
  LinTri *ltri;
  LinTet *ltet;
  LinLine *lline;
  LagrangeLine *lagLines[10];
  DoFManager DOFM;
  int dof_per_node;
  int n_els;
  
public:
  ElementFactory();
  ElementFactory(Mesh &M, const DoFManager &dofm);
  ~ElementFactory();

  Element & getElement(unsigned elnum); //return NULL if element not implemented

  const GenericMesh<MTYPE,NSD> & getMesh() {return M;}
  inline unsigned get_n_dof() { return dof_per_node*(M.n_nodes()); }
  inline unsigned get_dof_per_node() { return dof_per_node; }
};

YAFEL_NAMESPACE_CLOSE

#endif //#ifndef _YAFEL_ELEMENTFACTORY_HPP
