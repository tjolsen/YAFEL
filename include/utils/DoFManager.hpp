#ifndef _YAFEL_DOFMANAGER_HPP
#define _YAFEL_DOFMANAGER_HPP

#include "yafel_globals.hpp"
#include "yafel_typedefs.hpp"
#include "new_mesh/Mesh.hpp"
#include "element/ElementType.hpp"
#include <vector>


YAFEL_NAMESPACE_OPEN

class DoFManager {

public:
    enum class ManagerType {
        CG,
        DG
    };


    DoFManager(const Mesh &M, ManagerType m_type, int polyOrder, int dof_per_node=1);

    void getGlobalDofs(int elnum, std::vector<int> &container);

//private:
    ElementType CellType_to_ElementType(CellType ct, int polyOrder) const;
    void make_cg_dofs(const Mesh &M);
    void make_dg_dofs(const Mesh &M);
    int make_raw_dofs(const Mesh &M);
    coordinate<> interpolate_from_corners(coordinate<> xlocal,
                                          const std::vector<coordinate<>> &corners,
                                          CellType ct) const noexcept;
    void recombine_all_duplicates();


    int dof_per_node;
    int polyOrder;
    ManagerType managerType;
    std::vector<coordinate<>> dof_nodes;
    std::vector<int> element_offsets;
    std::vector<int> elements;

};



/*
class DoFManager {
  
public:
  using size_type = std::size_t;
  using element_container = std::vector<size_type>;

  DoFManager() : DoFManager(1) {}
  DoFManager(size_type dofpn) : _dof_per_node(dofpn) {}

  inline size_type dof_per_node() const {return _dof_per_node;}

  virtual size_type global_index(size_type elnum, size_type local_node, size_type comp) const = 0;
  virtual size_type n_dofs() const = 0;

protected:
  size_type _dof_per_node;
};
*/
YAFEL_NAMESPACE_CLOSE

#endif
