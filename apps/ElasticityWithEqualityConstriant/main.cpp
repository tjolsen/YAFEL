//
// Created by tyler on 2/7/18.
//

#include <lin_alg/linear_solvers/solvers/EigenBICGSTAB.hpp>
#include "yafel.hpp"
#include "LinearElasticity.hpp"


class KKTMatrix;

using Eigen::SparseMatrix;

/*
namespace Eigen {
namespace internal {

template<>
struct traits<MatrixReplacement> : public Eigen::internal::traits<Eigen::SparseMatrix<double, Eigen::RowMajor> >
{
};

}
}
*/

/**
 * A type extending the Eigen library to hold a block sparse matrix:
 *
 * KKT = [A B^T]
 *       [B  0 ]
 */
/*
class KKTMatrix : public Eigen::EigenBase<KKTMatrix>
{
public:
    // Required typedefs, constants, and method:
    typedef double Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = Eigen::SparseMatrix<double, Eigen::RowMajor>::IsRowMajor
    };

    auto rows() const { return A->rows() + B->rows(); }

    auto cols() const { return A->cols() + B->rows(); }

    //my api
    using InternalMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

    KKTMatrix(InternalMatrix const &a, InternalMatrix const &b) : A(&a), B(&b) {
        assert(a.rows() == a.cols());
        assert(a.cols() == b.cols());
    }

    auto& getA() const { return *A; }
    auto& getB() const { return *B; }

private:
    const InternalMatrix *A, *B;
};

namespace Eigen {
namespace internal {
template<typename Rhs>
struct generic_product_impl<KKTMatrix, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
        : generic_product_impl_base<KKTMatrix, Rhs, generic_product_impl<MatrixReplacement, Rhs> >
{
    typedef typename Product<MatrixReplacement, Rhs>::Scalar Scalar;

    template<typename Dest>
    static void scaleAndAddTo(Dest &dst, const MatrixReplacement &lhs, const Rhs &rhs, const Scalar &alpha)
    {
        // This method should implement "dst += alpha * lhs * rhs" inplace,
        // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
        assert(alpha == Scalar(1) && "scaling is not implemented");

        dst.topRows(lhs.getA().rows()).noalias() += lhs.getA()*rhs.topRows(lhs.getA().cols())
                                                    + lhs.getB().transpose()*rhs.bottomRows(lhs.getB().rows());

        dst.bottomRows(lhs.getB().rows()).noalias() += lhs.getB()*rhs.topRows(lhs.getA().cols());
    }
};
}
}
*/

using namespace yafel;


int main()
{
    constexpr int nsd = 3;
    Mesh M("mesh.msh");
    int p = 1;


    DoFManager dofm(M, DoFManager::ManagerType::CG, p, nsd);
    FESystem feSystem(dofm);

    CGAssembly<LinearElasticity<nsd>>(feSystem);

    //std::cout << "Assembly done\n";

    //Ordinary dirichlet BC on surface 1: u = (.01, 0, 0)
    DirichletBC bc0(dofm, 0.001, 0);
    bc0.selectByRegionID(1);
    DirichletBC bc1(dofm, 0, 1);
    bc1.selectByRegionID(1);
    DirichletBC bc2(dofm, 0, 2);
    bc2.selectByRegionID(1);

    DirichletBC bc3(dofm, 0, 2);
    bc2.selectByRegionID(3);

    auto &K = feSystem.getGlobalTangent();
    auto &rhs = feSystem.getGlobalResidual();

    //std::cout << "Dirichlet BCs done\n";

    //Constraint: u \cdot (1,1,0) == 0 on surface 2
    int ConstraintSurface = 2;
    std::vector<Eigen::Triplet<double>> constraint_triplets;
    std::vector<int> node_container;
    Tensor<nsd, 1> surfaceNormal{1, 1, 0};
    //surfaceNormal = surfaceNormal/norm(surfaceNormal);

    struct normal_constraint
    {
        int constrained_dof;
        std::vector<int> dofs;
        std::vector<double> coefficient;
        double inhomogeneity{0.0};
    };
    std::vector<normal_constraint> constraints;
    for (int c = 0; c < M.nCells(); ++c) {
        if (dofm.cell_region_idx[c] == ConstraintSurface) {
            dofm.getGlobalNodes(c, node_container);

            for (auto n : node_container) {
                normal_constraint C;
                C.constrained_dof = n * nsd;
                C.dofs = {(n*nsd)+1};//, n * nsd + 2};
                C.coefficient = {-surfaceNormal(1)};//, surfaceNormal(2)};
                constraints.push_back(C);
            }
        }
    }
    //std::cout << "Constraint creation done\n";


    //Condense linear system for constraints
    //rhs condensation:
    for (auto &C : constraints) {
        auto q_i = C.constrained_dof;
        for (int l = 0; l < C.dofs.size(); ++l) {
            rhs(C.dofs[l]) += C.coefficient[l] * rhs(q_i);
        }
        rhs(q_i) = 0;
    }
    std::vector<int> constraint_entry(rhs.rows(), -1);
    {
        int idx{0};
        for (auto &C : constraints) {
            constraint_entry[C.constrained_dof] = idx;
            ++idx;
        }
    }

    //In this instance, I know I'm not adding any nonzeros, so I can just loop over the current nonzeros
    for (int row = 0; row < rhs.rows(); ++row) {

        if (constraint_entry[row] == -1) {
            //UnConstrained DoF

            for (int idx = K.outerIndexPtr()[row]; idx < K.outerIndexPtr()[row + 1]; ++idx) {
                int col = K.innerIndexPtr()[idx];

                if (constraint_entry[col] != -1) {
                    //Constrained column in unconstrained row
                    auto entry = constraint_entry[col];
                    auto &value = K.valuePtr()[idx];
                    for (int l = 0; l < constraints[entry].dofs.size(); ++l) {
                        K.coeffRef(row, constraints[entry].dofs[l]) += value * constraints[entry].coefficient[l];
                    }

                    rhs(row) -= value * constraints[entry].inhomogeneity;
                    value = 0;
                }

            }

        } else {
            //Constrained dof, distribute row
            auto row_entry = constraint_entry[row];
            for (int idx = K.outerIndexPtr()[row]; idx < K.outerIndexPtr()[row + 1]; ++idx) {
                auto &value = K.valuePtr()[idx];
                auto col = K.innerIndexPtr()[idx];

                if (constraint_entry[col] == -1) {
                    //unconstrained column, set old entry to zero, after distributing values
                    for (int l = 0; l < constraints[row_entry].dofs.size(); ++l) {
                        K.coeffRef(constraints[row_entry].dofs[l], col) +=
                                value * constraints[row_entry].coefficient[l];
                    }
                    value = 0;
                } else {
                    //constrained column
                    //
                    auto col_entry = constraint_entry[col];
                    for (int l = 0; l < constraints[row_entry].dofs.size(); ++l) {
                        for (int q = 0; q < constraints[col_entry].dofs.size(); ++q) {

                            K.coeffRef(constraints[row_entry].dofs[l],
                                    constraints[col_entry].dofs[q])
                                    +=
                                    value
                                    *constraints[row_entry].coefficient[l]
                                    *constraints[col_entry].coefficient[q];

                        }

                        rhs(constraints[row_entry].dofs[l]) -=
                                value*constraints[row_entry].coefficient[l]*constraints[col_entry].inhomogeneity;
                    }
                    if(row == col)
                        value = 10000;
                    else
                        value = 0;
                }


            } //end loop over row nonzeros

            //fix up vector
            for(int l=0; l<constraints[row_entry].dofs.size(); ++l) {
                rhs(constraints[row_entry].dofs[l]) += rhs(row)*constraints[row_entry].coefficient[l];
            }
            rhs(row) = 0;

        } //end if row is constrained


    }

    //std::cout << "Condensation done\n";


    bc0.apply(K, rhs);
    //bc1.apply(K, rhs);
    bc2.apply(K, rhs);
    //bc3.apply(K, rhs);

    auto solverTag = LinearSolve::VCLConjugateGradientTag{};
    feSystem.getSolution() = LinearSolve::solve(K, rhs,solverTag);
    //std::cout << "Solve done\n";


    std::function<void(FESystem &, OutputFrame &)> captureFunc = [](FESystem &feSys, OutputFrame &frame) {
        frame.time = 0;
        OutputData::DataLocation dataLocation = OutputData::DataLocation::Point;
        std::string sol_name = "U";
        OutputData::DataType dt = OutputData::DataType::Vector;
        auto dat = std::make_shared<OutputData>(feSys.getSolution(),
                                                sol_name, dataLocation, dt, std::vector<int>(nsd, 1));
        frame.addData(dat);
    };

    SimulationOutput simulationOutput("output", BackendType::VTU);
    simulationOutput.captureFrame(feSystem, captureFunc);

}