//
// Created by tyler on 3/12/17.
//

#ifndef YAFEL_MESH_TYPEDEFS_HPP
#define YAFEL_MESH_TYPEDEFS_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/Tensor.hpp"

YAFEL_NAMESPACE_OPEN

/**
 * \class CellType
 * \brief Enum denoting geometric type of a cell in a mesh
 *
 * This enum class defines the types of cells that can be
 * represented in a mesh. Often, these will be linear tris, tets,
 * quads, hexes, and lines. If necessary, higher-order
 * versions will be added, if mesh generation of curved boundaries
 * is desired.
 *
 * In addition there will be (already is?) a set of functions
 * to convert from Gmsh element types (or other mesh generators)
 * to CellTypes, and from CellTypes to VTK Element Types.
 * This will be necessary when importing/exporting meshes
 * from external utilities.
 *
 * The naming convention is: <CellTopology><NumberOfNodes>.
 * More concretely, names will be of the form:
 *
 * - Line2
 * - Line3
 * - Tri3
 * - Tri6
 * - Quad4
 * - Quad9
 * - Hex8
 * - Hex20
 * - ...
 */
enum class CellType : int {
    Line2,
    Tri3,
    Tet4,
    Quad4,
    Hex8,
    Point1,  // Useful for MPM?
    None
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_MESH_TYPEDEFS_HPP
