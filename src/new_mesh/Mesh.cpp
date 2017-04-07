//
// Created by tyler on 3/15/17.
//

#include "new_mesh/Mesh.hpp"
#include <fstream>

YAFEL_NAMESPACE_OPEN


Mesh::Mesh(DefinitionScheme definitionScheme,
           const std::vector<coordinate<>> &geometryNodes,
           const std::vector<int> &cellNodes,
           const std::vector<int> &cellOffsets,
           const std::vector<CellType> &cellTypes,
           const std::vector<std::vector<int>> &cellTags)
        : definitionScheme_(definitionScheme),
          geometryNodes_(geometryNodes),
          cellNodes_(cellNodes),
          cellOffsets_(cellOffsets),
          cellTypes_(cellTypes),
          cellTags_(cellTags),
          internal_faces_(0) {}


Mesh::Mesh(const std::string &fname)
        : definitionScheme_(DefinitionScheme::Explicit)
{
    parse_gmsh(fname);
}

void Mesh::getCellNodes(int elem, std::vector<int> &container) const
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