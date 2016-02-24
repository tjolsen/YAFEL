#include "utils/DirBC.hpp"

YAFEL_NAMESPACE_OPEN

//========================================================================
//DirBC::DirBC(const GenericMesh &m, const DoFManager &dofm, std::size_t tagID,
//	     std::size_t comp, const SpatialFunction<double> &sfunc)
//{
//
//  this->comp = comp;
//  std::size_t nNodes = m.get_n_nodes();
//  std::size_t nElems = m.get_n_elems();
//  std::size_t ndofs = dofm.n_dofs(nNodes);
//  
//  bcmask.clear();
//  bcmask.resize(ndofs, false);
//  
//  for(std::size_t e=0; e<nElems; ++e) {
//    std::size_t id = m.el_tags[e][0];
//    if(id == tagID) {
//      for(std::size_t n=0; n<m.elements[e].size(); ++n) {
//	std::size_t node = m.elements[e][n];
//	std::size_t index = dofm.index(node, comp);
//	bcdofs.push_back(index);
//	bcmask[index] = true;
//	double val = sfunc(m.nodal_coords[node]);
//	bcvals.push_back(val);
//      }
//    }
//  }
//  
//}
//
////========================================================================
//DirBC::DirBC(const std::vector<std::size_t> bcdofs_,
//	     const std::vector<double> bcvals_,
//	     const std::vector<bool> bcmask_) :
//  bcdofs(bcdofs_), bcvals(bcvals_), bcmask(bcmask_)
//{}
//
//========================================================================
//void DirBC::apply(sparse_csr &Ksys, Vector<double> &Fsys) {
//
//  Vector Ubc(Fsys.size(), 0.0);
//  for(std::size_t i=0; i<bcvals.size(); ++i) {
//    Ubc(bcdofs[i]) = bcvals[i];
//  }
//
//  //store for later use
//  this->ubc = Ubc;
//
//  Fsys -= Ksys*Ubc;
//
//  for(std::size_t i=0; i<bcvals.size(); ++i) {
//    Fsys(bcdofs[i]) = bcvals[i];
//  }
//  
//  //set Ksys entries to 0 or 1 to apply BC's
//  for(std::size_t r=0; r<Ksys.getRows(); ++r) {
//    for(std::size_t i=Ksys.row_ptr[r]; i<Ksys.row_ptr[r+1]; ++i) {
//      std::size_t c = Ksys.col_index[i];
//      if(bcmask[r] || bcmask[c])
//	Ksys.data[i] = 0.0;
//      if(bcmask[r] && bcmask[c])
//	Ksys.data[i] = 1.0;
//    }
//  }
//  
//}
//
//========================================================================
//========================================================================
//========================================================================
YAFEL_NAMESPACE_CLOSE
