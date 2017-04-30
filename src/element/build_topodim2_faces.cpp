//
// Created by tyler on 4/23/17.
//

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>

YAFEL_NAMESPACE_OPEN

static std::vector<std::vector<int>> build_topodim2_faces(const std::vector<coordinate<>> &xi_all,
                                                          const std::function<bool(coordinate<>)> &isBoundary,
                                                          double theta0)
{
    using std::cos;
    using std::sin;
    double pi = std::atan(1.0) * 4;
    auto Theta = [pi](coordinate<> xi) { return std::atan2(xi(1), xi(0)); };

    Tensor<3,2> R;
    R(0,0) = cos(theta0);
    R(0,1) = -sin(theta0);
    R(1,0) = sin(theta0);
    R(1,1) = cos(theta0);
    R(2,2) = 1;
    coordinate<> xi0{0.25, 0.25, 0};

    std::vector<double> thetaVec;
    std::vector<int> bnd_idxs;
    int idx = 0;
    for (auto xi : xi_all) {
        if (isBoundary(xi)) {
            bnd_idxs.push_back(idx);
            double th = Theta(R*(xi-xi0));
            thetaVec.push_back(th);
        }
        ++idx;
    }

    std::vector<int> tmp_idx(bnd_idxs.size());
    idx = 0;
    for (auto &i : tmp_idx) {
        i = idx++;
    }

    auto sort_func_f = [&thetaVec](int l, int r) { return thetaVec[l] < thetaVec[r]; };
    auto sort_func_r = [&thetaVec](int l, int r) { return thetaVec[l] > thetaVec[r]; };

    std::vector<int> f_bnd(tmp_idx);
    std::sort(f_bnd.begin(), f_bnd.end(), sort_func_f);
    for (auto &f : f_bnd) {
        f = bnd_idxs[f];
    }

    idx = 0;
    std::vector<int> r_bnd(tmp_idx);
    std::sort(r_bnd.begin(), r_bnd.end(), sort_func_r);
    for (auto &r : r_bnd) {
        r = bnd_idxs[r];
    }

    return {f_bnd, r_bnd};
}


void Element::build_quad_faces()
{
    double pi = std::atan(1.0) * 4;

    face_perm.push_back({build_topodim2_faces(localMesh.getGeometryNodes(),
                                              [](auto xi) { return std::abs(xi(1) + 1) < 1.0e-6; }, pi / 2)});

    face_perm.push_back({build_topodim2_faces(localMesh.getGeometryNodes(),
                                              [](auto xi) { return std::abs(xi(0) - 1) < 1.0e-6; }, 0)});

    face_perm.push_back({build_topodim2_faces(localMesh.getGeometryNodes(),
                                              [](auto xi) { return std::abs(xi(1) - 1) < 1.0e-6; },
                                              -pi / 2)});

    face_perm.push_back({build_topodim2_faces(localMesh.getGeometryNodes(),
                                              [](auto xi) { return std::abs(xi(0) + 1) < 1.0e-6; }, pi)});
}


void Element::build_tri_faces()
{
    double pi = std::atan(1.0) * 4;
    face_perm.push_back({build_topodim2_faces(
            localMesh.getGeometryNodes(),
            [](auto xi) { return std::abs(xi(1)) < 1.0e-6; },
            -pi / 2
    )});

    face_perm.push_back({build_topodim2_faces(
            localMesh.getGeometryNodes(),
            [](auto xi) { return std::abs(1 - sum(xi)) < 1.0e-6; },
            pi / 4
    )});

    face_perm.push_back({build_topodim2_faces(
            localMesh.getGeometryNodes(),
            [](auto xi) { return std::abs(xi(0)) < 1.0e-6; },
            -pi

    )});

}
YAFEL_NAMESPACE_CLOSE