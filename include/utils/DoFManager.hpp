#ifndef __YAFEL_DOFMANAGER_HPP
#define __YAFEL_DOFMANAGER_HPP

#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

class DoFManager {
  
private:
  unsigned dofPerNode;

public:
  DoFManager() : DoFManager(1) {}
  DoFManager(unsigned dofpn) : dofPerNode(dofpn) {}
  inline unsigned index(unsigned node, unsigned comp) const {return node*dofPerNode + comp;}
  inline unsigned getDofPerNode() const {return dofPerNode;}
  inline unsigned n_dofs(unsigned nNodes) const {return dofPerNode*nNodes;}
};

YAFEL_NAMESPACE_CLOSE

#endif
