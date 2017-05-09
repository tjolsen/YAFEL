//
// Created by tyler on 5/2/17.
//

#include "element/ElementFactory.hpp"
#include "assembly/DGAssembly.hpp"
#include "fe_system/FESystem.hpp"
#include "output/SimulationOutput.hpp"
#include "time_integration/DGRK4.hpp"

using namespace yafel;

template<int NSD>
struct AdvectionPhysics
{
    std::vector<Eigen::PartialPivLU<Eigen::MatrixXd>> inverse_mass_matrices;
    bool mass_constructed{false};

    static constexpr int nsd() { return NSD; }

    static double initial_condition(const coordinate<> &x)
    {

        coordinate<> dx = x - coordinate<>{1.25, 1};

        int mask = dx(0)>0 && dx(0)<.5 && dx(1) > 0 && dx(1) <.5;
        double pi = 4*std::atan(1.0);
        return mask*std::sin(2*pi*dx(0))*std::sin(2*pi*dx(1));

    }

    static inline Tensor<NSD, 1> convection_velocity(coordinate<> x, double)
    {

        //double omega = 3.14159;
        //coordinate<> dx = x - coordinate<>{1, 1, x(2)};

        //double r = norm(dx);
        //double theta = std::atan2(dx(1), dx(0));

        //Tensor<NSD, 1> v{-r * omega * std::sin(theta), r * omega * std::cos(theta)};

        //return Tensor<NSD,1>{-dx(1)/r, dx(0)/r};

        //return v;
        return {-.75, .75};
    };

    static inline double source(const coordinate<> &, double)
    {
        return 0.0;
    }

    template<typename TU, typename TR>
    static inline void LocalResidual(const Element &E, int qpi,
                              coordinate<> xqp, double time,
                              Eigen::MatrixBase<TU> &U_el,
                              Eigen::MatrixBase<TR> &R_el)
    {
        auto c = convection_velocity(xqp, time);

        double U{0};
        for (int A = 0; A < U_el.rows(); ++A) {
            U += U_el(A) * E.shapeValues[qpi](A);
        }

        for (int i = 0; i < R_el.rows(); ++i) {
            auto gradW = make_TensorMap<NSD, 1>(&E.shapeGrad(i, 0));
            auto val = (U * dot(gradW, c)
                        + source(xqp, time) * E.shapeValues[qpi](i));
            R_el(i) += val * E.jxw;
        }
    }

    template<typename T>
    static void LocalMass(const Element &E, int qpi, double,
                          Eigen::MatrixBase<T> &M_el)
    {
        M_el += (E.jxw * E.shapeValues[qpi]) * E.shapeValues[qpi].transpose();
    }


    static double BoundaryFlux(Tensor<NSD, 1> n, coordinate<> x, double t, double U)
    {
        auto vdotn = dot(n, convection_velocity(x, t));
        if (vdotn > 0) {
            return vdotn * U;
        } else {
            return vdotn * U;
        }
    }

    static inline double Flux(Tensor<NSD, 1> n, coordinate<> x, double t, double Uplus, double Uminus)
    {
        auto vdotn = dot(n, convection_velocity(x, t));
        return 0.5 * (vdotn * (Uplus + Uminus) + std::abs(vdotn) * (Uplus - Uminus));
    }

};


int main()
{
    constexpr int NSD = 2;
    Mesh M("mesh.msh");
    M.buildInternalFaces();
    int p = 3;
    int dofpn = 1;

    AdvectionPhysics<NSD> AP;

    double Tfinal = 1;
    double dt = .001;
    int output_freq = 1;


    DoFManager dofm(M, DoFManager::ManagerType::DG, p, dofpn);

    FESystem feSystem(dofm);

    auto &U = feSystem.getSolution();

    //set initial condition
    for (int i = 0; i < dofm.dof_nodes.size(); ++i) {
        U(i) = AdvectionPhysics<NSD>::initial_condition(dofm.dof_nodes[i]);
    }

    SimulationOutput simulationOutput("output", BackendType::HDF5);

    simulationOutput.captureFrame(feSystem);

    DGRK4 dgrk4(dt);

    int NT = static_cast<int>(Tfinal / dt);
    for (int ti = 1; ti <= NT; ++ti) {
        std::cout << ti << " / " << NT << std::endl;
        dgrk4.step(feSystem,AP);

        if (ti % output_freq == 0) {
            simulationOutput.captureFrame(feSystem);
        }
    }

}