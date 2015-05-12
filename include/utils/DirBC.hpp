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

public:
  unsigned comp;
  std::vector<unsigned> bcdofs;
  std::vector<double> bcvals;
  std::vector<bool> bcmask;
  DoFManager DOFM;
  Vector ubc;


  DirBC(const Mesh &m, const DoFManager &dofm, unsigned tagID, 
	       unsigned comp, const SpatialFunction<double> &sfunc);

  DirBC(const std::vector<unsigned> bcdofs_,
	       const std::vector<double> bcvals_,
	       const std::vector<bool> bcmask_);
  
  void apply(sparse_csr &Ksys, Vector &Fsys);
  inline Vector getUbc() const {return ubc;}
};

YAFEL_NAMESPACE_CLOSE

#endif
