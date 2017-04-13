cmake_minimum_required(VERSION 3.5)
project(YAFEL)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math -Wall -Wextra -DNDEBUG -march=native")

set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -fno-omit-frame-pointer -fsanitize=address")

set(SOURCE_FILES
        include/element/Element.hpp
        include/element/ElementFactory.hpp
        include/element/ElementType.hpp
        include/element/ShapeFunctionUtils.hpp

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
        include/lin_alg/tensor/tensor_functions/spectral_decomposition.hpp
        include/lin_alg/tensor/tensor_functions/tensor_dot.hpp
        include/lin_alg/tensor/tensor_functions/UnaryFunctions.hpp
        include/lin_alg/tensor/tensor_functions/update_assignment.hpp


        include/mesh/ElementIterator.hpp
        include/mesh/mesh_typedefs.hpp
        include/mesh/Mesh.hpp
        include/mesh/CellFace.hpp

        include/output/OutputBackend.hpp
        include/output/OutputData.hpp
        include/output/OutputFrame.hpp
        include/output/OutputMesh.hpp
        include/output/VTUBackend.hpp
        #include/output/MatrixVisualization.hpp
        #include/output/VTKDGMesh.hpp
        #include/output/VTKMesh.hpp
        #include/output/VTKObject.hpp
        #include/output/VTKOutput.hpp
        #include/output/VTKScalarData.hpp
        #include/output/VTKTensorData.hpp
        #include/output/VTKTimeOutput.hpp
        #include/output/VTKVectorData.hpp

        include/quadrature/LobattoPoints1D.hpp
        include/quadrature/QuadratureRule.hpp

        include/utils/ConcurrentQueue.hpp
        include/utils/old_DirBC.hpp
        include/utils/DoFManager.hpp
        include/utils/DualNumber.hpp
        include/utils/ElementVtkType.hpp
        include/utils/Range.hpp
        include/utils/ScalarTraits.hpp
        include/utils/SpatialFunction.hpp
        include/utils/TaskScheduler.hpp
        include/utils/ThreadPool.hpp
        include/yafel.hpp
        include/yafel_globals.hpp
        include/yafel_typedefs.hpp


        src/element/Element.cpp
        src/element/ElementFactory.cpp
        src/element/ShapeFunctionUtils.cpp
        src/quadrature/LobattoPoints1D.cpp
        src/mesh/Mesh.cpp
        src/mesh/build_faces.cpp
        src/mesh/parse_gmsh.cpp

        src/output/OutputBackend.cpp
        src/output/OutputData.cpp
        src/output/VTUBackend.cpp

        #src/output/VTKObject.cpp
        #src/output/VTKOutput.cpp
        #src/output/VTKScalarData.cpp
        #src/output/VTKTimeOutput.cpp

        src/quadrature/QuadratureRule.cpp
        src/quadrature/TriangleQuadratureRules.cpp

        src/utils/DoFManager.cpp
        #src/utils/old_DirBC.cpp
        #src/utils/ElementVtkType.cpp
        )

include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(yafel ${SOURCE_FILES})

add_executable(sandbox experimental/sandbox.cpp)
target_link_libraries(sandbox yafel)


# Add any subdirectories to the project
add_subdirectory(apps)