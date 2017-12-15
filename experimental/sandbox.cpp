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

        K_el += E.shapeGrad*E.shapeGrad.transpose();

    }


};


int main()
{

    double L = 1;

    Mesh M(Mesh::DefinitionScheme::Explicit);
    std::vector<coordinate<>> nodes{{0,0,0}, {L,0,0}, {0,L,0}, {0,0,L}};
    std::vector<int> cells{0,1,2,3};
    std::vector<int> offsets{0,4};
    std::vector<CellType> ctypes{CellType::Tet4};

    M.setGeometryNodes(nodes);
    M.setCellNodes(cells);
    M.setOffsets(offsets);
    M.setCellTypes(ctypes);
    M.setCellTags({{0}});

    int p = 1;
    int dofpn = 1;
    DoFManager dofm(M, DoFManager::ManagerType::CG,p,dofpn);

    constexpr int NSD = 3;

    ElementFactory EF;
    auto et = dofm.element_types[0];

    auto & E = EF.getElement(et);
    E.update<NSD>(0,0,dofm);

    cout << E.shapeGrad << endl << endl;
    cout << E.shapeValues[0] << endl << endl;
    cout << E.jxw << endl << endl;

    FESystem system(dofm, NSD);

    CGAssembly<Laplacian<NSD>>(system);

    auto &R = system.getGlobalResidual();
    double s = 0;
    for(int i=0; i<R.rows(); ++i) {
        s += R(i);
    }

    cout << system.getGlobalResidual() << endl << endl << "s = " << s << endl;
    cout << system.getGlobalTangent() << endl << endl;
    return 0;
}


