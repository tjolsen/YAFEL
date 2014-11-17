#ifndef _YAFEL_ELEMENTFACTORY_HPP
#define _YAFEL_ELEMENTFACTORY_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "element/LinQuad.hpp"
#include "element/LinTri.hpp"
#include "element/LinTet.hpp"
#include "mesh/Mesh.hpp"

YAFEL_NAMESPACE_OPEN

class ElementFactory {

private:
  Mesh *Mp;
  LinQuad *lq;
  LinTri *ltri;
  LinTet *ltet;
  int dof_per_node;
  int n_els;
  
public:
  ElementFactory();
  ElementFactory(Mesh &M, int dofPerNode);
  ~ElementFactory();

  Element *getElement(int elnum); //return NULL if element not implemented

  Mesh *getMesh() {return Mp;}
  inline int get_n_dof() { return dof_per_node*(Mp->get_n_nodes()); }
  inline int get_dof_per_node() { return dof_per_node; }
};

YAFEL_NAMESPACE_CLOSE

#endif //#ifndef _YAFEL_ELEMENTFACTORY_HPP
