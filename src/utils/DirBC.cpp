#include "utils/DirBC.hpp"

YAFEL_NAMESPACE_OPEN

DirBC::DirBC(const std::vector<unsigned> &bcn, const std::vector<unsigned> &bcc,
	     const std::vector<double> &bcv, const DoFManager &dofm) 
  : bcnodes(bcn), bccomps(bcc), bcvals(bcv), DOFM(dofm)
{}

void DirBC::apply(sparse_csr &Ksys, Vector &Fsys) {
  
  
  
}


YAFEL_NAMESPACE_CLOSE
