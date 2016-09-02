#ifndef _YAFEL_DIRBC_HPP
#define _YAFEL_DIRBC_HPP

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
#include "lin_alg/operators.hpp"
#include "utils/DoFManager.hpp"
#include "utils/SpatialFunction.hpp"
#include <vector>
#include <iostream>

YAFEL_NAMESPACE_OPEN

class DirBC {

public:
    using size_type = std::size_t;

    std::vector<size_type> bcdofs;
    std::vector<double> bcvals;
    std::vector<bool> bcmask;
    Vector<double> ubc;
  
    template<unsigned NSD>
    DirBC(const GmshMesh<NSD> &m, const DoFManager &dofm, size_type tagID, 
	  size_type comp, const SpatialFunction<NSD,double> &sfunc);

    DirBC(const std::vector<size_type> bcdofs_,
	  const std::vector<double> bcvals_,
	  const std::vector<bool> bcmask_)
	: bcdofs(bcdofs_), bcvals(bcvals_), bcmask(bcmask_), ubc(bcmask_.size(),0)
	{};
  
    void apply(sparse_csr<double> &Ksys, Vector<double> &Fsys);
    inline Vector<double> get_ubc() const {return ubc;}
};



/*
 * Implementation starts here
 */

template<unsigned NSD>
DirBC::DirBC(const GmshMesh<NSD> &M, const DoFManager &dofm, size_type tagID, 
             size_type comp, const SpatialFunction<NSD,double> &sfunc) 
    : bcdofs(), 
      bcvals(), 
      bcmask(dofm.n_dofs(),false), 
      ubc(dofm.n_dofs(),0)
{

    //size_type nNodes = M.get_n_nodes();
    size_type nElems = M.n_elements();
    //size_type ndofs = dofm.n_dofs();
  
    for(size_type e=0; e<nElems; ++e) {
	size_type id = M._element_tags[e][0];
	if(id == tagID) {
	    for(size_type n=0; n<M.element(e).size(); ++n) {
		size_type index = dofm.global_index(e, n, comp);
		bcdofs.push_back(index);
		bcmask[index] = true;
		size_type nodenum = M.element(e)[n];
		double val = sfunc(M.node(nodenum));
		bcvals.push_back(val);
	    }
	}
    }
  
}


YAFEL_NAMESPACE_CLOSE

#endif
