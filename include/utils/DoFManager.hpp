#ifndef _YAFEL_DOFMANAGER_HPP
#define _YAFEL_DOFMANAGER_HPP

#include "yafel_globals.hpp"
#include "yafel_typedefs.hpp"
#include "mesh/Mesh.hpp"
#include "element/ElementType.hpp"
#include <vector>


YAFEL_NAMESPACE_OPEN

class DoFManager
{

public:
    //Denote DoF layout: CG, DG, (HDG, EDG?), (HDG_Skeleton, EDG_Skeleton?), (Composite?)
    enum class ManagerType
    {
        CG,
        DG
    };


    DoFManager(const Mesh &M, ManagerType m_type, int polyOrder, int dof_per_node = 1);

    void getGlobalDofs(int elnum, std::vector<int> &container) const;

    void getGlobalNodes(int elnum, std::vector<int> &container) const;

    //void getGlobalFaceDofs(int elnum, std::vector<int> &container) const;
    //void getGlobalFaceNodes(int elnum, std::vector<int> &container) const;

    ElementType CellType_to_ElementType(CellType ct, int polyOrder) const;

    inline int nCells() const { return static_cast<int>(element_types.size()); }

    int dof_per_node;
    int polyOrder;
    ManagerType managerType;
    std::vector<coordinate<>> dof_nodes;
    std::vector<int> element_offsets;
    std::vector<int> elements;
    std::vector<ElementType> element_types;
    std::vector<int> cell_region_idx;

    std::vector<CellFace> interior_faces;
    std::vector<int> face_offsets;
    std::vector<int> face_left_local_nodes;
    std::vector<int> face_right_local_nodes;


    std::vector<int> mesh_corner_idxs;
    std::vector<int> mesh_corner_offsets;
private:
    void make_cg_dofs(const Mesh &M);

    void make_dg_dofs(const Mesh &M);

    int make_raw_dofs(const Mesh &M, bool match_face_dofs = false);

    coordinate<> interpolate_from_corners(coordinate<> xlocal,
                                          const std::vector<coordinate<>> &corners,
                                          CellType ct) const noexcept;

    void recombine_all_duplicates();

};

YAFEL_NAMESPACE_CLOSE

#endif
