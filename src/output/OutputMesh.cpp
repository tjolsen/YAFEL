//
// Created by tyler on 5/1/17.
//

#include "output/OutputMesh.hpp"
#include "element/ElementFactory.hpp"

YAFEL_NAMESPACE_OPEN

OutputMesh::OutputMesh(const DoFManager &dofm_)
        : dofm(&dofm_)
{
    //Need an element factory to build visualization topology
    ElementFactory EF(1);

    int n_points = dofm->dof_nodes.size();
    int n_parent_cells = dofm->element_offsets.size() - 1;

    local_cells_per_cell.resize(n_parent_cells);

    int offset{0};
    std::vector<int> local_cell;
    std::vector<int> global_nodes;
    for (auto e : IRange(0, n_parent_cells)) {
        auto et = dofm->element_types[e];
        if (et.elementTopology == ElementTopology::None) {
            continue;
        }
        auto &E = EF.getElement(et);
        dofm->getGlobalNodes(e, global_nodes);

        local_cells_per_cell[e] = E.localMesh.nCells();

        for (auto lc : IRange(0, E.localMesh.nCells())) {
            auto et = dofm->CellType_to_ElementType(E.localMesh.getCellType(lc), 1);
            E.localMesh.getCellNodes(lc, local_cell);
            expanded_cell_offsets.push_back(offset);
            offset += local_cell.size();

            expanded_cell_element_type.push_back(et);

            for (auto n : local_cell) {
                expanded_cells.push_back(global_nodes[n]);
            }
        }
    }
    expanded_cell_offsets.push_back(offset);
}


YAFEL_NAMESPACE_CLOSE