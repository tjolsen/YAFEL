#ifndef _YAFEL_CG_DOFMANAGER_HPP
#define _YAFEL_CG_DOFMANAGER_HPP

/*
 * CG_DoFManager:
 *
 * Class to handle the mapping of element number, local node number,
 * and node DoF component to global DoF numbering. The primary purpose
 * for this is during the assembly process, when local element residual
 * vectors and tangent matrices are assembled into the global structures.
 *
 * The class is designed for use in continuous Galerkin FEM codes.
 * As such, each node holds dof_per_node degrees of freedom, and the
 * global degree of freedom can be computed directly from a global node
 * number (taken from the mesh) without the need for additional data storage.
 */

#include "yafel_globals.hpp"
#include "mesh/GenericMesh.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

template<typename MT, unsigned NSD>
class CG_DoFManager : public DoFManager {

public:
  using size_type = typename DoFManager::size_type;

  CG_DoFManager(const GenericMesh<MT,NSD> &m, size_type dofpn) : DoFManager(dofpn), M(m) {}

  size_type global_index(size_type elnum, size_type local_node, size_type component=0) const {
    return M.element(elnum)[local_node]*_dof_per_node + component;
  }

  size_type n_dofs() const {
    return M.n_nodes()*_dof_per_node;
  }

private:
  const GenericMesh<MT,NSD> & M;
};

YAFEL_NAMESPACE_CLOSE

#endif
