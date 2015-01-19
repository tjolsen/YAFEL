#ifndef _DIRBC_HPP
#define _DIRBC_HPP

#include "yafel_globals.hpp"
#include "mesh/Mesh.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"
#include "utils/SpatialFunction.hpp"
#include <vector>


YAFEL_NAMESPACE_OPEN

class DirBC {

private:
  unsigned comp;
  std::vector<unsigned> bcnodes;
  std::vector<double> bcvals;
  std::vector<bool> bcmask;
  DoFManager DOFM;

public:
  DirBC::DirBC(const Mesh &m, const DoFManager &dofm, unsigned tagID, 
	       unsigned comp, const SpatialFunction<double> &sfunc);
  
  void apply(sparse_csr &Ksys, Vector &Fsys);
};

YAFEL_NAMESPACE_CLOSE

#endif
