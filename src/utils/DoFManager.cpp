//
// Created by tyler on 4/7/17.
//

#include "utils/DoFManager.hpp"
#include "utils/Range.hpp"
#include "element/ElementFactory.hpp"

YAFEL_NAMESPACE_OPEN

DoFManager::DoFManager(const Mesh &M,
                       ManagerType m_type,
                       int polyOrder,
                       int dof_per_node)
        : dof_per_node(dof_per_node),
          polyOrder(polyOrder),
          managerType(m_type),
          mesh_corner_idxs(M.getCellVector()),
          mesh_corner_offsets(M.getOffsetVector())
{

    cell_region_idx.resize(M.nCells());
    for (auto c : IRange(0, M.nCells())) {
        cell_region_idx[c] = M.getCellTags(c)[0];
    }

    switch (m_type) {
        case ManagerType::CG:
            make_cg_dofs(M);
            break;
        case ManagerType::DG:
            interior_faces = M.getInternalFaces();
            make_dg_dofs(M);
            break;
    }
}


void DoFManager::getGlobalDofs(int elnum, std::vector<int> &container) const
{
    int n_dofs = dof_per_node * (element_offsets[elnum + 1] - element_offsets[elnum]);

    container.clear();
    container.reserve(n_dofs);

    for (auto idx : IRange(element_offsets[elnum], element_offsets[elnum + 1])) {
        int node = elements[idx];
        for (auto dof : IRange(0, dof_per_node)) {
            container.push_back(node * dof_per_node + dof);
        }
    }
}

void DoFManager::getGlobalNodes(int elnum, std::vector<int> &container) const
{
    int n_dofs = element_offsets[elnum + 1] - element_offsets[elnum];

    container.clear();
    container.reserve(n_dofs);

    for (auto idx : IRange(element_offsets[elnum], element_offsets[elnum + 1])) {
        int node = elements[idx];
        container.push_back(node);
    }
}

void DoFManager::getGlobalFaceNodes(int fnum, std::vector<int> &container) const
{
    int ndofs = face_offsets[fnum+1] - face_offsets[fnum];
    container.clear();
    container.reserve(ndofs);

    for(auto idx : IRange(face_offsets[fnum], face_offsets[fnum+1])) {
        container.push_back(idx);
    }
}

void DoFManager::getGlobalFaceDofs(int fnum, std::vector<int> &container) const
{

    int ndofs = dof_per_node * (face_offsets[fnum+1] - face_offsets[fnum]);
    container.clear();
    container.reserve(ndofs);

    for(auto idx : IRange(face_offsets[fnum], face_offsets[fnum+1])) {
        for(auto dof : IRange(0,dof_per_node)) {
            container.push_back(dof_per_node*idx + dof);
        }
    }
}

void DoFManager::getLocalFaceNodes(
        int fnum,
        std::vector<int> &container,
        const std::vector<int> &source_container) const
{
    int ndofs = face_offsets[fnum+1] - face_offsets[fnum];
    container.clear();
    container.reserve(ndofs);

    for(auto idx : IRange(face_offsets[fnum], face_offsets[fnum+1])) {
        container.push_back(source_container[idx]);
    }
}


ElementType DoFManager::CellType_to_ElementType(CellType ct, int polyOrder) const
{

    ElementTopology topology;
    int topodim{0};

    switch (ct) {
        case CellType::Line2:
            topodim = 1;
            topology = ElementTopology::TensorProduct;
            break;
        case CellType::Tri3:
            topodim = 2;
            topology = ElementTopology::Simplex;
            break;
        case CellType::Tet4:
            topodim = 3;
            topology = ElementTopology::Simplex;
            break;
        case CellType::Quad4:
            topodim = 2;
            topology = ElementTopology::TensorProduct;
            break;
        case CellType::Hex8:
            topodim = 3;
            topology = ElementTopology::TensorProduct;
            break;
        default:
            topology = ElementTopology::None;
            break;
    }

    return ElementType(topology, topodim, polyOrder);
}


void DoFManager::make_cg_dofs(const Mesh &M)
{
    if(polyOrder > 1) {
        //make a raw mesh
        make_raw_dofs(M);

        //recombine duplicate nodes and update triangulation
        recombine_all_duplicates();
    }

    else if (polyOrder == 1) {
        this->elements = M.getCellVector();
        this->element_offsets = M.getOffsetVector();
        this->dof_nodes = M.getGeometryNodes();
        this->element_types.resize(M.nCells());
        for(int i=0; i<M.nCells(); ++i) {
            element_types[i] = CellType_to_ElementType(M.getCellType(i), polyOrder);
        }
    } else {
        throw std::runtime_error("DoFManager: Invalid polynomial order");
    }
}

void DoFManager::make_dg_dofs(const Mesh &M)
{
    //make a raw mesh
    make_raw_dofs(M);

    //Match nodes on mesh faces
    match_face_nodes();
}


int DoFManager::make_raw_dofs(const Mesh &M)
{

    int ncells = M.nCells();
    int8_t max_td{0};

    ElementFactory EF(1);
    int offset{0};

    std::vector<coordinate<>> corner_coords(8);
    std::vector<int> corner_idxs;
    auto const &geom_nodes = M.getGeometryNodes();

    element_offsets.resize(ncells + 1);
    element_types.resize(ncells);

    int next_node{0};
    for (auto c : IRange(0, ncells)) {

        auto et = CellType_to_ElementType(M.getCellType(c), polyOrder);
        element_types[c] = et;
        if (et.elementTopology == ElementTopology::None) {
            continue;
        }

        max_td = std::max(max_td, et.topoDim);


        // corner nodes
        M.getCellNodes(c, corner_idxs);
        corner_coords.clear();
        corner_coords.reserve(corner_idxs.size());

        for (auto i : corner_idxs) {
            corner_coords.push_back(geom_nodes[i]);
        }

        auto &E = EF.getElement(et);

        int nLocalNodes = E.localMesh.nNodes();

        for (auto n : E.localMesh.getGeometryNodes()) {
            dof_nodes.push_back(interpolate_from_corners(n, corner_coords, M.getCellType(c)));
            elements.push_back(next_node++);
        }


        element_offsets[c] = offset;
        offset += nLocalNodes;
    }
    element_offsets[ncells] = offset;

    return max_td;
}

void DoFManager::match_face_nodes()
{

    int nFaces = interior_faces.size();
    std::vector<int> left_global_nodes;
    std::vector<int> right_global_nodes;

    ElementFactory EF(1);


    for (int f = 0; f < nFaces; ++f) {
        face_offsets.push_back(face_left_local_nodes.size());
        auto &F = interior_faces[f];
        int left_elem = F.left;
        int left_flocal = F.left_flocal;

        int right_elem = F.right;
        int right_flocal = F.right_flocal;

        if (left_elem < 0) {
            auto &ER = EF.getElement(element_types[right_elem]);
            for (auto nr : ER.face_nodes[right_flocal]) {
                face_left_local_nodes.push_back(-1);
                face_right_local_nodes.push_back(nr);
            }
            continue;

        } else if (right_elem < 0) {
            auto &EL = EF.getElement(element_types[left_elem]);
            for (auto nl : EL.face_nodes[left_flocal]) {
                face_left_local_nodes.push_back(nl);
                face_right_local_nodes.push_back(-1);
            }
            continue;
        }


        getGlobalNodes(left_elem, left_global_nodes);
        getGlobalNodes(right_elem, right_global_nodes);

        auto &EL = EF.getElement(element_types[left_elem]);
        auto &ER = EF.getElement(element_types[right_elem]);

        auto &left_local_nodes = EL.face_nodes[left_flocal];
        auto &right_local_nodes = ER.face_nodes[right_flocal];


        for (auto nl : left_local_nodes) {
            face_left_local_nodes.push_back(nl);
            auto &xl = dof_nodes[left_global_nodes[nl]];

            for (auto nr : right_local_nodes) {
                auto &xr = dof_nodes[right_global_nodes[nr]];
                if (norm(xl - xr) < 1.0e-6) {
                    face_right_local_nodes.push_back(nr);
                    break;
                }
            }
        }

    }
    face_offsets.push_back(face_left_local_nodes.size());

}

coordinate<> DoFManager::interpolate_from_corners(coordinate<> xlocal,
                                                  const std::vector<coordinate<>> &corners,
                                                  CellType ct) const noexcept
{

    switch (ct) {
        case CellType::Tri3:
            return (1 - xlocal(0) - xlocal(1)) * corners[0]
                   + xlocal(0) * corners[1]
                   + xlocal(1) * corners[2];

        case CellType::Tet4:
            return (1 - xlocal(0) - xlocal(1) - xlocal(2)) * corners[0]
                   + xlocal(0) * corners[1]
                   + xlocal(1) * corners[2]
                   + xlocal(2) * corners[3];

        case CellType::Quad4:
            return (1.0 / 4.0) * (
                    (1 - xlocal(0)) * (1 - xlocal(1)) * corners[0] +
                    (1 + xlocal(0)) * (1 - xlocal(1)) * corners[1] +
                    (1 + xlocal(0)) * (1 + xlocal(1)) * corners[2] +
                    (1 - xlocal(0)) * (1 + xlocal(1)) * corners[3]
            );

        case CellType::Hex8:
            return (1.0 / 8.0) * (
                    (1 - xlocal(0)) * (1 - xlocal(1)) * (1 - xlocal(2)) * corners[0] +
                    (1 + xlocal(0)) * (1 - xlocal(1)) * (1 - xlocal(2)) * corners[1] +
                    (1 + xlocal(0)) * (1 + xlocal(1)) * (1 - xlocal(2)) * corners[2] +
                    (1 - xlocal(0)) * (1 + xlocal(1)) * (1 - xlocal(2)) * corners[3] +
                    (1 - xlocal(0)) * (1 - xlocal(1)) * (1 + xlocal(2)) * corners[4] +
                    (1 + xlocal(0)) * (1 - xlocal(1)) * (1 + xlocal(2)) * corners[5] +
                    (1 + xlocal(0)) * (1 + xlocal(1)) * (1 + xlocal(2)) * corners[6] +
                    (1 - xlocal(0)) * (1 + xlocal(1)) * (1 + xlocal(2)) * corners[7]
            );


        default:
            return coordinate<>{};
    }

}

void DoFManager::recombine_all_duplicates()
{

    int N = dof_nodes.size();
    std::vector<int> idxs(N);
    coordinate<> boundingBox_min{dof_nodes[0]};
    coordinate<> boundingBox_max{dof_nodes[0]};
    for (auto i : IRange(0, N)) {
        idxs[i] = i;
        boundingBox_min = min(boundingBox_min, dof_nodes[i]);
        boundingBox_max = max(boundingBox_max, dof_nodes[i]);
    }

    double max_dim = max(boundingBox_max - boundingBox_min);
    double snap_value = max_dim * 1.0e-6;

    auto snap = [snap_value](auto x) -> coordinate<> { return round(x / snap_value) * snap_value; };


    auto sort_func = [this, snap](int L, int R) {
        const auto xL = snap(dof_nodes[L]);
        const auto xR = snap(dof_nodes[R]);
        return std::lexicographical_compare(
                xL.data.rbegin(), xL.data.rend(),
                xR.data.rbegin(), xR.data.rend()
        );
    };
    //auto uniq_func_idx = [this](int L, int R) { return norm(this->dof_nodes[L] - this->dof_nodes[R]) < 1.0e-6; };
    auto uniq_func_x = [](auto x, auto y) { return norm(x - y) < 1.0e-6; };

    //lexicographical sort of nodes
    std::sort(idxs.begin(), idxs.end(), sort_func);

    std::vector<coordinate<>> sorted_nodes(N);
    std::vector<int> reverse_idxs(N);
    for (auto i : IRange(0, N)) {
        sorted_nodes[i] = dof_nodes[idxs[i]];
        reverse_idxs[idxs[i]] = i;
    }

    std::vector<int> uniq_idx(N);
    std::vector<int> uniq_uniq_idx(1, 0);
    uniq_idx[0] = 0;
    for (auto i : IRange(1, N)) {
        if (uniq_func_x(sorted_nodes[i], sorted_nodes[i - 1])) {
            uniq_idx[i] = uniq_idx[i - 1];
        } else {
            uniq_uniq_idx.push_back(i);
            uniq_idx[i] = i;
        }
    }

    std::vector<int> reverse_u_idx(uniq_idx.size());
    int current{0};
    reverse_u_idx[0] = current;
    for (auto i : IRange(1, static_cast<int>(uniq_idx.size()))) {
        if (uniq_idx[i] == uniq_idx[i - 1]) {
            reverse_u_idx[i] = current;
        } else {
            reverse_u_idx[i] = ++current;
        }
    }


    std::vector<coordinate<>> unique_nodes(uniq_uniq_idx.size());
    for (auto i : IRange(0, static_cast<int>(uniq_uniq_idx.size()))) {
        unique_nodes[i] = sorted_nodes[uniq_uniq_idx[i]];
    }

    for (auto &n : elements) {
        n = reverse_u_idx[uniq_idx[reverse_idxs[n]]];
    }

    dof_nodes.swap(unique_nodes);
}

YAFEL_NAMESPACE_CLOSE