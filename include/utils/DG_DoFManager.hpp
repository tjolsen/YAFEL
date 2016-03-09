#ifndef __YAFEL_DG_DOFMANAGER_HPP
#define __YAFEL_DG_DOFMANAGER_HPP


/*
 * DG_DoFManager
 *
 * Class to handle mapping of element number, local node number,
 * and node DoF component to a global DoF numbering. This is used
 * during assembly of global vectors and matrices.
 *
 * This class is designed for use with Discontinuous Galerkin
 * codes, and as such will assume that each element has
 * nodes_per_element*dof_per_node degrees of freedom. It
 * performs a pre-processing step by running through the mesh
 * and finding the starting global DoF for each element.
 * This is necessary to handle heterogeneous (tri/quad or tet/hex)
 * meshes, and to handle gmsh-style meshes where boundary elements
 * are included in the mesh, but not integrated.
 */

#include "yafel_globals.hpp"
#include "mesh/GenericMesh.hpp"
#include "utils/DoFManager.hpp"
#include <vector>
#include <exception>

YAFEL_NAMESPACE_OPEN

template<typename MT, unsigned NSD>
class DG_DoFManager : public DoFManager {

public:
  using size_type = typename DoFManager::size_type;
  
  DG_DoFManager(const GenericMesh<MT,NSD> &m, size_type dofpn);

  size_type global_index(size_type elnum, size_type local_node, size_type component=0) const {
    if(component >= _dof_per_node) {
      throw std::domain_error("Component exceeds dof_per_node");
    }
    if(local_node >= M.element(elnum).size()) {
      throw std::domain_error("local_node exceeds number of nodes in element");
    }
    
    return dof_offset[elnum] + local_node*_dof_per_node + component;
  }

  size_type n_dofs() const {
    return dof_offset[M.n_elements()];
  }
  
private:
  const GenericMesh<MT,NSD> &M;
  std::vector<size_type> dof_offset;
};



//==========================================================
/*
 * Implementation
 */
//==========================================================

/*
 * dof_offsets is the first global dof associated with an element.
 * Only the elements whose topological dimension is equal to the spatial dimension
 * of the DG_DoFManager.
 */
template<typename MT, unsigned NSD>
DG_DoFManager<MT,NSD>::DG_DoFManager(const GenericMesh<MT,NSD> &m, size_type dofpn)
  : DoFManager(dofpn), 
    M(m), 
    dof_offset(m.n_elements()+1,size_type(0))
{
  size_type offset=0;
  for(size_type e=0; e<m.n_elements(); ++e) {
    
    dof_offset[e] = offset;
    
    switch(m.element_type(e)) {
    case ElementType::LINEAR_LINE:
    case ElementType::QUADRATIC_LINE:
    case ElementType::CUBIC_LINE:
      if(NSD==1) {
        offset += dofpn*m.element(e).size();
        break;
      }
    case ElementType::LINEAR_TRI:
    case ElementType::QUADRATIC_TRI:
    case ElementType::CUBIC_TRI:
    case ElementType::LINEAR_QUAD:
    case ElementType::QUADRATIC_QUAD:
    case ElementType::CUBIC_QUAD:
      if(NSD==2) {
        offset += dofpn*m.element(e).size();
        break;
      }
    case ElementType::LINEAR_TET:
    case ElementType::QUADRATIC_TET:
    case ElementType::CUBIC_TET:
    case ElementType::LINEAR_HEX:
    case ElementType::QUADRATIC_HEX:
    case ElementType::CUBIC_HEX:
      if(NSD==3) {
        offset += dofpn*m.element(e).size();
        break;
      }
    case ElementType::DG_QUAD:
    case ElementType::NULL_ELEMENT:
      break;
    }//end switch
    
    dof_offset[m.n_elements()] = offset;

  }//end for
  
}

YAFEL_NAMESPACE_CLOSE

#endif
