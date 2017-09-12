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


// Assembly routines
#include "assembly/CGAssembly.hpp"
#include "assembly/DGAssembly.hpp"
#include "assembly/LocalSmoothingGradient.hpp"
#include "assembly/ZZGradientRecovery.hpp"

#include "boundary_conditions/DirichletBC.hpp"
#include "fe_system/FESystem.hpp"

#include "quadrature/LobattoPoints1D.hpp"
#include "quadrature/QuadratureRule.hpp"

// Utils
#include "utils/DoFManager.hpp"
#include "utils/DualNumber.hpp"
#include "utils/parallel/TaskScheduler.hpp"
#include "utils/parallel/parfor.hpp"
#include "utils/parallel/wait_all.hpp"
#include "utils/parallel/wait_any.hpp"

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
#include "output/SimulationOutput.hpp"
//#include "output/MatrixVisualization.hpp"

#endif
