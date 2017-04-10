//
// Created by tyler on 3/11/17.
//

#ifndef YAFEL_TENSORS_HPP
#define YAFEL_TENSORS_HPP
/**
 * \file
 * \brief Include the tensor library
 *
 * Note, this header does not include the metaprogramming headers.
 * These are included as needed in individual files in order to
 * keep people from inadvertently poking around and changing something.
 */

#include "yafel_globals.hpp"

//Base TensorExpression type. Contains necessary forward declarations
#include "lin_alg/tensor/TensorExpression.hpp"

//Expression types
#include "lin_alg/tensor/tensor_expression_types/BinaryOperations.hpp"
#include "lin_alg/tensor/tensor_expression_types/UnaryOperations.hpp"
#include "lin_alg/tensor/tensor_expression_types/Tensor.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorContraction.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorCwiseBinaryOp.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorCwiseUnaryOp.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorDyadicProduct.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorFunctor.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorMap.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorScaled.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorSlice.hpp"


//Functions and operators
#include "lin_alg/tensor/tensor_functions/BinaryFunctions.hpp"
#include "lin_alg/tensor/tensor_functions/UnaryFunctions.hpp"
#include "lin_alg/tensor/tensor_functions/operators.hpp"
#include "lin_alg/tensor/tensor_functions/tensor_dot.hpp"
#include "lin_alg/tensor/tensor_functions/update_assignment.hpp"
#include "lin_alg/tensor/tensor_functions/rank2_specializations.hpp"
#include "lin_alg/tensor/tensor_functions/spectral_decomposition.hpp"


#endif //YAFEL_TENSORS_HPP
