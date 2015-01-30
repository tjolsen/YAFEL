#include "utils/DirBC.hpp"

YAFEL_NAMESPACE_OPEN

//========================================================================
DirBC::DirBC(const Mesh &m, const DoFManager &dofm, unsigned tagID,
	     unsigned comp, const SpatialFunction<double> &sfunc)
{

  this->comp = comp;
  unsigned nNodes = m.get_n_nodes();
  unsigned nElems = m.get_n_elems();
  unsigned ndofs = dofm.getNDofs(nNodes);
  
  bcmask.clear();
  bcmask.resize(ndofs, false);
  
  for(unsigned e=0; e<nElems; ++e) {
    unsigned id = m.el_tags[e][0];
    if(id == tagID) {
      for(unsigned n=0; n<m.elements[e].size(); ++n) {
	unsigned node = m.elements[e][n];
	unsigned index = dofm.index(node, comp);
	bcdofs.push_back(index);
	bcmask[index] = true;
	double val = sfunc(m.nodal_coords[node]);
	bcvals.push_back(val);
      }
    }
  }
  
}

//========================================================================
void DirBC::apply(sparse_csr &Ksys, Vector &Fsys) {
  
  Vector Ubc(Fsys.getLength(), 0.0);
  for(unsigned i=0; i<bcvals.size(); ++i) {
    Ubc(bcdofs[i]) = bcvals[i];
  }

  //store for later use
  this->ubc = Ubc;

  Fsys -= Ksys*Ubc;

  for(unsigned i=0; i<bcvals.size(); ++i) {
    Fsys(bcdofs[i]) = bcvals[i];
  }
  
  //set Ksys entries to 0 or 1 to apply BC's
  for(unsigned r=0; r<Ksys.getRows(); ++r) {
    for(unsigned i=Ksys.row_ptr[r]; i<Ksys.row_ptr[r+1]; ++i) {
      unsigned c = Ksys.col_index[i];
      if(bcmask[r] || bcmask[c])
	Ksys.data[i] = 0.0;
      if(bcmask[r] && bcmask[c])
	Ksys.data[i] = 1.0;
    }
  }
  
}

//========================================================================
//========================================================================
//========================================================================
YAFEL_NAMESPACE_CLOSE
