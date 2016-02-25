#ifndef __YAFEL_DOFMANAGER_HPP
#define __YAFEL_DOFMANAGER_HPP

#include "yafel_globals.hpp"
#include <cstddef>

YAFEL_NAMESPACE_OPEN

class DoFManager {
  
public:
  using size_type = std::size_t;

  DoFManager() : DoFManager(1) {}
  DoFManager(size_type dofpn) : _dof_per_node(dofpn) {}

  inline size_type dof_per_node() const {return _dof_per_node;}

  virtual size_type global_index(size_type elnum, size_type local_node, size_type comp) const = 0;
  virtual size_type n_dofs() const = 0;

protected:
  size_type _dof_per_node;
};

YAFEL_NAMESPACE_CLOSE

#endif
