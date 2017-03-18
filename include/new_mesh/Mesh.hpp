//
// Created by tyler on 3/12/17.
//

#ifndef YAFEL_MESH_HPP
#define YAFEL_MESH_HPP

#include "yafel_globals.hpp"
#include "yafel_typedefs.hpp"
#include "mesh_typedefs.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

/**
 * \class Mesh
 * \brief Interface for meshes.
 *
 * This is the base class for meshes.
 * Specialized mesh types (Structured, unstructured, rectilinear...)
 * may inherit from this to provide specific functionality or storage savings.
 *
 * Geometry may be specified implicitly or explicitly,
 * and this can be used to affect the behavior of the
 * member functions. (eg: choosing whether to explicitly store
 * geometric/topological information or to generate it on the fly)
 *
 */
class Mesh
{
public:
    enum class DefinitionScheme : int
    {
        Explicit,
        Implicit
    };

    /**
     * \brief Construct a mesh with a specified storage scheme
     * @param definitionScheme
     * @param geometryNodes node coordinates
     * @param cellNodes nodes comprising cells
     * @param cellOffsets offset of each cell in cellNodes vector
     * @param cellTypes Types of each cell
     */
    Mesh(DefinitionScheme definitionScheme,
         std::vector<coordinate<>> geometryNodes = {},
         std::vector<int> cellNodes = {},
         std::vector<int> cellOffsets = {},
         std::vector<CellType> cellTypes = {});


    /**
     * \brief Build internal mesh faces. Implementation and storage TBD.
     */
    void buildInternalFaces();

    /**
     * Get node indices of a cell
     * @param cellnum Cell index
     * @param container vector to fill with node indices
     *
     * "container" is modified to hold cell nodes. No reallocation
     * is done unless the container needs to grow. Will be resized
     * to allow easy iteration over nodes for caller.
     */
    void getCellNodes(int cellnum, std::vector<int> &container) const;


    /**
     * Setters for internal structures
     */
    inline void setGeometryNodes(std::vector<coordinate<>> &&gn)
    { geometryNodes_ = gn; }

    inline void setCellNodes(std::vector<int> &&cn)
    { cellNodes_ = cn; }

    inline void setOffsets(std::vector<int> &&off)
    { cellOffsets_ = off; }

    inline void setCellTypes(std::vector<CellType> &&ct)
    { cellTypes_ = ct; }


    /**
     * Getters for internal structures
     */
    const std::vector<coordinate<>> &getGeometryNodes() const
    { return geometryNodes_; }


    //Get number of cells
    inline virtual int nCells() const noexcept
    {
        return cellTypes_.size();
    }

protected:
    DefinitionScheme definitionScheme_;

    // Containers for explicit geometry node storage
    // Stored in a "compressed" format, with element
    // "elem" comprising nodes from [elementOffsets_[elem], elementOffsets_[elem+1])
    // (note half-open interval).
    // cellOffsets_ has length "num_elems+1", with the last element
    // cellOffsets_[num_elems] = num_nodes. (Same as sparse matrix compressed storage)
    std::vector<coordinate<>> geometryNodes_;
    std::vector<int> cellNodes_;
    std::vector<int> cellOffsets_;
    std::vector<CellType> cellTypes_;

    virtual inline void getCellNodesImplicit(int cell, std::vector<int> &container) const noexcept
    {}
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_MESH_HPP
