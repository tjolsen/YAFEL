cmake_minimum_required(VERSION 3.5)
project(YAFEL)

find_package(OpenMP REQUIRED)
find_package(OpenCL REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)
find_package(HDF5 COMPONENTS C HL REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})


set(OPT_FLAGS
        "-march=native -malign-data=cacheline -ffast-math -fopenmp -finline-limit=1000000")

set(VIENNACL_CXX_FLAGS "-DVIENNACL_WITH_OPENCL")
#set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VIENNACL_CXX_FLAGS} -O3 -Wall -Wno-sign-compare -DNDEBUG -DEIGEN_DEFAULT_INDEX_TYPE=int ${OPT_FLAGS}")

set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -fno-omit-frame-pointer")

set(SOURCE_FILES
        include/assembly/AssemblyRequirement.hpp
        include/assembly/CGAssembly.hpp
        include/assembly/DGAssembly.hpp
        include/assembly/LocalSmoothingGradient.hpp
        include/assembly/ZZGradientRecovery.hpp

        include/boundary_conditions/DirichletBC.hpp

        include/element/Element.hpp
        include/element/ElementFactory.hpp
        include/element/ElementType.hpp
        include/element/ShapeFunctionUtils.hpp
        include/element/element_boundary_nodes.hpp

        include/fe_system/FESystem.hpp

        include/lin_alg/tensor/tensors.hpp
        include/lin_alg/tensor/TensorExpression.hpp
        include/lin_alg/tensor/mp_utils/contraction_mp_utils.hpp
        include/lin_alg/tensor/mp_utils/map_mp_utils.hpp
        include/lin_alg/tensor/mp_utils/mp_tuple_utils.hpp
        include/lin_alg/tensor/mp_utils/sequence_functions.hpp
        include/lin_alg/tensor/mp_utils/sequences.hpp
        include/lin_alg/tensor/mp_utils/slice_mp_utils.hpp
        include/lin_alg/tensor/mp_utils/TypeList.hpp
        include/lin_alg/tensor/tensor_expression_types/BinaryOperations.hpp
        include/lin_alg/tensor/tensor_expression_types/Tensor.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorContraction.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorCwiseBinaryOp.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorCwiseUnaryOp.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorDyadicProduct.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorFunctor.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorMap.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorPermutation.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorScaled.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorSlice.hpp
        include/lin_alg/tensor/tensor_expression_types/UnaryOperations.hpp
        include/lin_alg/tensor/tensor_functions/BinaryFunctions.hpp
        include/lin_alg/tensor/tensor_functions/operators.hpp
        include/lin_alg/tensor/tensor_functions/rank2_specializations.hpp
        include/lin_alg/tensor/tensor_functions/rank1_specializations.hpp
        include/lin_alg/tensor/tensor_functions/Reductions.hpp
        include/lin_alg/tensor/tensor_functions/spectral_decomposition.hpp
        include/lin_alg/tensor/tensor_functions/tensor_dot.hpp
        include/lin_alg/tensor/tensor_functions/UnaryFunctions.hpp
        include/lin_alg/tensor/tensor_functions/update_assignment.hpp


        include/mesh/ElementIterator.hpp
        include/mesh/mesh_typedefs.hpp
        include/mesh/Mesh.hpp
        include/mesh/CellFace.hpp

        include/output/HDFBackend.hpp
        include/output/OutputBackend.hpp
        include/output/OutputData.hpp
        include/output/OutputFrame.hpp
        include/output/OutputMesh.hpp
        include/output/SimulationOutput.hpp
        include/output/VTUBackend.hpp

        include/quadrature/LobattoPoints1D.hpp
        include/quadrature/QuadratureRule.hpp

        include/time_integration/DGRK4.hpp

        include/utils/BasicTimer.hpp
        include/utils/ConcurrentQueue.hpp
        include/utils/DoFManager.hpp
        include/utils/DualNumber.hpp
        include/utils/ElementVtkType.hpp
        include/utils/Range.hpp
        include/utils/ScalarTraits.hpp
        include/utils/SmallVector.hpp
        include/utils/SpatialFunction.hpp
        include/utils/parallel/TaskScheduler.hpp
        include/utils/parallel/Task.hpp
        include/utils/parallel/wait_all.hpp
        include/utils/parallel/wait_any.hpp
        include/utils/parallel/parfor.hpp
        include/utils/parallel/ReductionVariable.hpp
        include/utils/ThreadPool.hpp

        include/yafel.hpp
        include/yafel_globals.hpp
        include/yafel_typedefs.hpp
        include/yafel_config.hpp
        include/utils/parallel/yafel_parallel.hpp

        src/boundary_conditions/DirichletBC.cpp

        src/element/Element.cpp
        src/element/ElementFactory.cpp
        src/element/ShapeFunctionUtils.cpp
        src/element/element_boundary.cpp
        #src/element/build_tet_faces.cpp
        #src/element/build_topodim2_faces.cpp
        src/element/make_tensorproduct_element.cpp
        src/element/make_simplex_element.cpp

        src/mesh/CellFace.cpp
        src/mesh/Mesh.cpp
        src/mesh/build_faces.cpp
        src/mesh/parse_gmsh.cpp

        src/output/HDFBackend.cpp
        src/output/OutputBackend.cpp
        src/output/OutputData.cpp
        src/output/OutputMesh.cpp
        src/output/SimulationOutput.cpp
        src/output/VTUBackend.cpp

        src/quadrature/LobattoPoints1D.cpp
        src/quadrature/QuadratureRule.cpp
        src/quadrature/TriangleQuadratureRules.cpp
        src/quadrature/TetrahedronQuadratureRules.cpp

        src/time_integration/DGRK4.cpp

        src/utils/DoFManager.cpp
        src/utils/parallel/TaskScheduler.cpp
        )

link_libraries("OpenCL")
link_libraries(${HDF5_C_LIBRARIES})
link_libraries(${HDF5_C_HL_LIBRARIES})

include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(yafel ${SOURCE_FILES})

add_executable(sandbox experimental/sandbox.cpp)
target_link_libraries(sandbox yafel)


# Add any subdirectories to the project
add_subdirectory(apps)