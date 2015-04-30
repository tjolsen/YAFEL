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

  for(unsigned i=0; i<10; ++i) {
    lagLines[i] = nullptr;
  }
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
    
    // Lagrange Line elements of increasing orders. Orders 7-10 not documented...
    // orders 11+ not supported by gmsh
  case 8: 
    if(lagLines[1]==nullptr)
      lagLines[1] = new LagrangeLine(DOFM,2);
    return lagLines[1];
  case 26:
    if(lagLines[2]==nullptr)
      lagLines[2] = new LagrangeLine(DOFM,3);
    return lagLines[2];
  case 27:
    if(lagLines[3]==nullptr)
      lagLines[3] = new LagrangeLine(DOFM,4);
    return lagLines[3];
  case 28:
    if(lagLines[4]==nullptr)
      lagLines[4] = new LagrangeLine(DOFM,5);
    return lagLines[4];
  case 62:
    if(lagLines[5]==nullptr)
      lagLines[5] = new LagrangeLine(DOFM,6);
    return lagLines[5];
  case 63:
    if(lagLines[6]==nullptr)
      lagLines[6] = new LagrangeLine(DOFM,7);
    return lagLines[6];
  case 64:
    if(lagLines[7]==nullptr)
      lagLines[7] = new LagrangeLine(DOFM,8);
    return lagLines[7];
  case 65:
    if(lagLines[8]==nullptr)
      lagLines[8] = new LagrangeLine(DOFM,9);
    return lagLines[8];
  case 66:
    if(lagLines[9]==nullptr)
      lagLines[9] = new LagrangeLine(DOFM,10);
    return lagLines[9];
  }

  return nullptr;
}




YAFEL_NAMESPACE_CLOSE
