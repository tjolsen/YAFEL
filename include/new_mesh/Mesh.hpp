//
// Created by tyler on 3/12/17.
//

#ifndef YAFEL_MESH_HPP
#define YAFEL_MESH_HPP

#include "yafel_globals.hpp"
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
     */
    Mesh(DefinitionScheme definitionScheme = DefinitionScheme::Explicit);


    /**
     * \brief Build internal mesh faces
     */
    void buildInternalFaces();


    inline void getElementNodes(int elem, std::vector<int> &container) const noexcept
    {
        if (definitionScheme_ == DefinitionScheme::Explicit) {
            int n_elem_nodes = elementOffsets_[elem + 1] - elementOffsets_[elem];
            container.clear();
            container.reserve(n_elem_nodes);

            for (int i = elementOffsets_[elem]; i < elementOffsets_[elem + 1]; ++i) {
                conatiner.push_back(elementNodes[i]);
            }
            return;
        } else {
            getElementNodesImplicit(elem, container);
            return;
        }
    }

protected:
    DefinitionScheme definitionScheme_;

    // Containers for explicit geometry node storage
    // Stored in a "compressed" format, with element
    // "elem" comprising nodes from [elementOffsets_[elem], elementOffsets_[elem+1])
    // (note half-open interval).
    // elementOffsets_ has length "num_elems+1", with the last element
    // elementOffsets_[num_elems] = num_nodes. (Same as sparse matrix compressed storage)
    std::vector<coordinate<>> geometryNodes_;
    std::vector<int> elementNodes;
    std::vector<int> elementOffsets_;


    virtual inline void getElementNodesImplicit(int element, std::vector<int> &container) {}
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_MESH_HPP
