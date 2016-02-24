#ifndef __YAFEL_DIRBC_HPP
#define __YAFEL_DIRBC_HPP

/*
 * Class to represent and facilitate the application of Dirichlet boundary conditions
 */


#include "yafel_globals.hpp"
#include "mesh/GenericMesh.hpp"
#include "mesh/GmshMesh.hpp"
#include "mesh/RectilinearMesh.hpp"
#include "lin_alg/access_sparse_matrix.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"
#include "utils/SpatialFunction.hpp"
#include <vector>


YAFEL_NAMESPACE_OPEN

class DirBC {

public:
  using size_type = std::size_t;

  std::vector<size_type> bcdofs;
  std::vector<double> bcvals;
  std::vector<bool> bcmask;
  DoFManager DOFM;
  Vector<double> ubc;
  
  template<unsigned NSD>
  DirBC(const GmshMesh<NSD> &m, const DoFManager &dofm, size_type tagID, 
        size_type comp, const SpatialFunction<NSD,double> &sfunc);

  DirBC(const std::vector<size_type> bcdofs_,
	       const std::vector<double> bcvals_,
	       const std::vector<bool> bcmask_)
    : bcdofs(bcdofs_), bcvals(bcvals_), bcmask(bcmask_)
  {};
  
  template<typename T>
  void apply(access_sparse_matrix<T,double> &Ksys, Vector<double> &Fsys);
  inline Vector<double> getUbc() const {return ubc;}
};



/*
 * Implementation starts here
 */

template<unsigned NSD>
DirBC(const GmshMesh<MT,NSD> &m, const DoFManager &dofm, size_type tagID, 
      size_type comp, const SpatialFunction<NSD,double> &sfunc) 
  : bcdofs(), 
    bcvals(), 
    bcmask(dofm.n_dofs(M.n_nodes()),false), 
    DOFM(dofm), 
    ubc(dofm.n_dofs(M.n_nodes()),0)
{

  size_type nNodes = m.get_n_nodes();
  size_type nElems = m.get_n_elems();
  size_type ndofs = dofm.n_dofs(nNodes);
  
  for(size_type e=0; e<nElems; ++e) {
    size_type id = m.el_tags[e][0];
    if(id == tagID) {
      for(size_type n=0; n<m.element(e).size(); ++n) {
	size_type nodenum = m.element(e)[n];
	size_type index = dofm.index(nodenum, comp);
	bcdofs.push_back(index);
	bcmask[index] = true;
	double val = sfunc(m.node(nodenum));
        bcvals.push_back(val);
      }
    }
  }
  
}


template<>
void apply(sparse_csr<double> &Kcsr, Vector<double> &Fsys) {

  for(std::size_t i=0; i<bcvals.size(); ++i) {
    ubc(bcdofs[i]) = bcvals[i];
  }

  Fsys -= Ksys*ubc;

  for(std::size_t i=0; i<bcvals.size(); ++i) {
    Fsys(bcdofs[i]) = bcvals[i];
  }
  
  //set Ksys entries to 0 or 1 to apply BC's
  for(std::size_t r=0; r<Ksys.rows(); ++r) {
    for(std::size_t i=Ksys.row_ptr[r]; i<Ksys.row_ptr[r+1]; ++i) {
      std::size_t c = Ksys.col_index[i];
      if(bcmask[r] || bcmask[c])
	Ksys._data[i] = 0.0;
      if(bcmask[r] && bcmask[c])
	Ksys._data[i] = 1.0;
    }
  }

}

YAFEL_NAMESPACE_CLOSE

#endif
