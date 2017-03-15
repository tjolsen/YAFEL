//
// Created by tyler on 3/15/17.
//

#include "new_mesh/Mesh.hpp"

YAFEL_NAMESPACE_OPEN


Mesh::Mesh(DefinitionScheme definitionScheme,
           std::vector<coordinate<>> geometryNodes,
           std::vector<int> cellNodes,
           std::vector<int> cellOffsets,
           std::vector<CellType> cellTypes)
        : definitionScheme_(std::move(definitionScheme)),
          geometryNodes_(std::move(geometryNodes)),
          cellNodes_(std::move(cellNodes)),
          cellOffsets_(std::move(cellOffsets)),
          cellTypes_(std::move(cellTypes))
{}


void Mesh::getCellNodes(int elem, std::vector<int> &container) const noexcept
{
    if (definitionScheme_ == DefinitionScheme::Explicit) {
        int n_cell_nodes = cellOffsets_[elem + 1] - cellOffsets_[elem];
        container.clear();
        container.reserve(n_cell_nodes);

        for (int i = cellOffsets_[elem]; i < cellOffsets_[elem + 1]; ++i) {
            container.push_back(cellNodes_[i]);
        }
        return;
    } else {
        getCellNodesImplicit(elem, container);
        return;
    }
}


YAFEL_NAMESPACE_CLOSE