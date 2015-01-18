#include "element/ElementFactory.hpp"

YAFEL_NAMESPACE_OPEN

ElementFactory::ElementFactory() {}

ElementFactory::ElementFactory(Mesh &M, const DoFManager &dofm) :
  Mp(&M), DOFM(dofm), dof_per_node(dofm.getDofPerNode())
{
  
  //initialize implementations of elements here... 
  lq = new LinQuad(DOFM);
  ltri = new LinTri(DOFM);
  ltet = new LinTet(DOFM);
  lline = new LinLine(DOFM);
}

ElementFactory::~ElementFactory() {
  delete lq;
  delete ltri;
  delete ltet;
  delete lline;
}

Element *ElementFactory::getElement(unsigned elnum) {
  
  int eltype = Mp->element_type[elnum];
  
  switch(eltype) {
  case 1: return lline; //2-node line
  case 2: return ltri; //3-node triangle
  case 3: return lq; //4-node quad
  case 4: return ltet; // 4-node tetrahedron (will be implemented)
  }

  return NULL;
}




YAFEL_NAMESPACE_CLOSE
