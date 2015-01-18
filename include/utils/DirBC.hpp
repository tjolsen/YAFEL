#ifndef _DIRBC_HPP
#define _DIRBC_HPP

#include "yafel_globals.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"
#include <vector>


YAFEL_NAMESPACE_OPEN

class DirBC {

private:
  std::vector<unsigned> bcnodes;
  std::vector<unsigned> bccomps;
  std::vector<double> bcvals;
  DoFManager DOFM;

public:
  DirBC(const std::vector<unsigned> &bcn, const std::vector<unsigned> &bcc,
	const std::vector<double> &bcv, const DoFManager &dofm);
  
  void apply(sparse_csr &Ksys, Vector &Fsys);
};

YAFEL_NAMESPACE_CLOSE

#endif
