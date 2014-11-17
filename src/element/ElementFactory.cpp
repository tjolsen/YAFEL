#include "element/ElementFactory.hpp"

YAFEL_NAMESPACE_OPEN

ElementFactory::ElementFactory() {}

ElementFactory::ElementFactory(Mesh &M, int dofPerNode) :
  Mp(&M), dof_per_node(dofPerNode)
{
  
  //initialize implementations of elements here... 
  lq = new LinQuad(dof_per_node);
  ltri = new LinTri(dof_per_node);
  ltet = new LinTet(dof_per_node);
}

ElementFactory::~ElementFactory() {
  delete lq;
  delete ltri;
}

Element *ElementFactory::getElement(int elnum) {
  
  int eltype = Mp->element_type[elnum];
  
  switch(eltype) {
  case 1: break; //2-node line
  case 2: return ltri; //3-node triangle
  case 3: return lq; //4-node quad
  case 4: return ltet; // 4-node tetrahedron (will be implemented)
  }

  return NULL;
}




YAFEL_NAMESPACE_CLOSE
