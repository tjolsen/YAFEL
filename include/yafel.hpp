#ifndef _YAFEL_HPP
#define _YAFEL_HPP

/*
  This file serves as the simplest interface to YAFEL by including everything.
  If desired, this file may be ignored in favor of more fine-grained control
  over the parts of the library that are included.
  It deliberately includes headers that are not meant to be user-facing
  as a convenience to library developers (aka me!).
 */

// Library-wide global declarations
#include "yafel_globals.hpp"

// Linear Algebra structures
#include "lin_alg/tensor/tensors.hpp" // <-- already a convenience header


// Utils
#include "utils/old_DirBC.hpp"
#include "utils/DoFManager.hpp"
#include "utils/DualNumber.hpp"
#include "utils/ElementType.hpp"
#include "quadrature/QuadratureRule.hpp"
#include "utils/SpatialFunction.hpp"

// Mesh
#include "mesh/mesh_typedefs.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/CellFace.hpp"

// Elements
#include "element/Element.hpp"
#include "element/ElementFactory.hpp"
#include "element/ElementType.hpp"
#include "element/ShapeFunctionUtils.hpp"

// Output
#include "output/VTKObject.hpp"
#include "output/VTKOutput.hpp"
#include "output/VTKScalarData.hpp"
#include "output/VTKVectorData.hpp"
#include "output/VTKTensorData.hpp"
#include "output/VTKMesh.hpp"
#include "output/VTKDGMesh.hpp"
#include "output/VTKTimeOutput.hpp"
#include "output/MatrixVisualization.hpp"

#endif
