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
          managerType(m_type)
{

    switch (m_type) {
        case ManagerType::CG:
            make_cg_dofs(M);
            break;
        case ManagerType::DG:
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
    //make a raw mesh
    make_raw_dofs(M);

    //recombine duplicate nodes and update triangulation
    recombine_all_duplicates();

}

void DoFManager::make_dg_dofs(const Mesh &M)
{
    //make a raw mesh
    make_raw_dofs(M);

}


int DoFManager::make_raw_dofs(const Mesh &M)
{

    int ncells = M.nCells();
    int8_t max_td{0};

    ElementFactory EF(1);
    int offset{0};

    std::vector<coordinate<>> corner_coords(8);
    std::vector<int> corner_idxs;
    auto const& geom_nodes = M.getGeometryNodes();

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
        std::cout << c << ": ";
        for(auto i : corner_idxs) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
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


coordinate<> DoFManager::interpolate_from_corners(coordinate<> xlocal,
                                                  const std::vector<coordinate<>> &corners,
                                                  CellType ct) const noexcept
{

    switch (ct) {
        case CellType::Tri3:
            return (1 - xlocal(0) - xlocal(1)) * corners[0] + xlocal(0) * corners[1] + xlocal(1) * corners[2];

        case CellType::Tet4:
            return (1 - xlocal(0) - xlocal(1) - xlocal(2)) * corners[0]
                   + xlocal(0) * corners[1] + xlocal(1) * corners[2]
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
            return coordinate<>();
    }

}


void DoFManager::recombine_all_duplicates()
{
    std::vector<int> idxs(dof_nodes.size(), 0);
    int idx = 0;
    for (auto &i : idxs)
        i = idx++;

    auto sort_func = [this](int L, int R) { return this->dof_nodes[L].data < this->dof_nodes[R].data; };
    auto uniq_func = [this](int L, int R) { return norm(this->dof_nodes[L] - this->dof_nodes[R]) < 1.0e-10; };


    //lexicographical sort of nodes
    std::sort(idxs.begin(), idxs.end(), sort_func);
    std::vector<int> reverse_idx(idxs.size());
    for(auto i : IRange(0,static_cast<int>(idxs.size()))) {
        reverse_idx[idxs[i]] = i;
    }

    std::vector<int> uniq_idx(idxs.size());
    uniq_idx[0] = idxs[0];
    for(auto i : IRange(1,static_cast<int>(idxs.size()))) {

        if(uniq_func(idxs[i], idxs[i-1])) {
            uniq_idx[i] = uniq_idx[i-1];
        }
        else {
            uniq_idx[i] = idxs[i];
        }
    }

    std::vector<int> unique_uniq_idx(uniq_idx);
    auto uuend = std::unique(unique_uniq_idx.begin(), unique_uniq_idx.end());
    unique_uniq_idx.resize(std::distance(unique_uniq_idx.begin(), uuend));
    std::vector<int> reverse_map(uniq_idx.size(),0);
    for(auto i : IRange(0,static_cast<int>(unique_uniq_idx.size()))) {
        reverse_map[unique_uniq_idx[i]] = i;
    }

    for(auto &n : elements) {
        n = reverse_map[uniq_idx[reverse_idx[n]]];
    }

    std::vector<coordinate<>> dof_nodes_unique;
    dof_nodes_unique.reserve(unique_uniq_idx.size());

    for(auto i : unique_uniq_idx) {
        dof_nodes_unique.push_back(dof_nodes[i]);
    }

    dof_nodes.swap(dof_nodes_unique);
}

YAFEL_NAMESPACE_CLOSE