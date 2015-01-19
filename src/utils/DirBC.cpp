#include "utils/DirBC.hpp"

YAFEL_NAMESPACE_OPEN

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
    unsigned id = m.elements[e].el_tags[0];
    if(id == tagID) {
      for(unsigned n=0; n<m.elements[e].size(); ++n) {
	unsigned node = m.elements[e][n];
	unsigned index = dofm.index(node, comp);
	bcnodes.push_back[node];
	bcmask[index] = true;
	double val = sfunc(m.nodal_coords[node]);
	bcvals.push_back(val);
      }
    }
  }
  
}

void DirBC::apply(sparse_csr &Ksys, Vector &Fsys) {
 
  
  
}


YAFEL_NAMESPACE_CLOSE
