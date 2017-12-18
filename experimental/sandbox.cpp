//
// Created by tyler on 3/14/17.
//

#include "yafel.hpp"
#include "PDE/PDEBase.hpp"
using namespace yafel;
using std::cout;
using std::endl;


template<int NSD>
struct Laplacian : PDEBase<NSD> {

    template<typename Utype, typename Rtype>
    static void LocalResidual(const Element& E, int qpi, coordinate<>, double,
                              Eigen::DenseBase<Utype> &u_el,
                              Eigen::DenseBase<Rtype> &R_el) {

        R_el += 6*E.shapeValues[qpi]*E.jxw;

    }

    template<typename Utype, typename Ktype>
    static void LocalTangent(const Element& E, int qpi, coordinate<>, double,
                              Eigen::DenseBase<Utype> &u_el,
                              Eigen::DenseBase<Ktype> &K_el) {

        K_el += E.shapeGrad*E.shapeGrad.transpose()*E.jxw;

    }


};


int main()
{

    double L = 1.0;

    Mesh M("mintetcube.msh");

    /*
    Mesh M(Mesh::DefinitionScheme::Explicit);
    std::vector<coordinate<>> nodes{{0,0,0}, {L,0,0}, {0,L,0}, {0,0,L}, {L,L,L}};
    std::vector<int> cells{0,1,2,3,1,2,3,4};
    std::vector<int> offsets{0,4,8};
    std::vector<CellType> ctypes{CellType::Tet4, CellType::Tet4};

    M.setGeometryNodes(nodes);
    M.setCellNodes(cells);
    M.setOffsets(offsets);
    M.setCellTypes(ctypes);
    M.setCellTags({{0}, {0}});
    */

    int p = 1;
    int dofpn = 1;
    DoFManager dofm(M, DoFManager::ManagerType::CG,p,dofpn);

    constexpr int NSD = 3;

    ElementFactory EF;
    auto et = dofm.element_types[0];

    auto & E = EF.getElement(et);
    E.update<NSD>(0,0,dofm);

    cout << E.shapeGradXi[0] << endl << endl;
    cout << E.shapeValues[0] << endl << endl;
    cout << E.jxw << endl << endl;

    FESystem system(dofm, NSD);
    CGAssembly<Laplacian<NSD>>(system);

    DirichletBC bc1(dofm,0.0);
    bc1.selectByFunction([&](auto &x){return std::abs(x(0)) < 1.0e-6; });
    //DirichletBC bc2(dofm,1.0);
    //bc2.selectByFunction([&](auto &x){return std::abs(x(0)-1) < 1.0e-6; });
    bc1.apply(system.getGlobalTangent(),system.getGlobalResidual());
    //bc2.apply(system.getGlobalTangent(),system.getGlobalResidual());

    std::cout << system.getGlobalTangent() << endl << endl;
    cout << system.getGlobalResidual() << endl << endl;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper | Eigen::Lower> solver;
    solver.compute(system.getGlobalTangent());
    Eigen::VectorXd result = solver.solve(system.getGlobalResidual());

    cout << "result = " << endl << result << endl  << endl;

    return 0;
}


