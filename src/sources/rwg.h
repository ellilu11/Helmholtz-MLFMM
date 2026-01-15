#pragma once

#include "source.h"
#include "triangle.h"

class RWG final : public Source {
public:
    friend class BC;

    RWG(std::shared_ptr<Excitation::PlaneWave>,
        size_t,
        const Eigen::Vector4i&, 
        const TriVec&);

    RWG(std::shared_ptr<Triangle>,
        std::shared_ptr<Triangle>);

    double getLeng() const { return leng; }

    vec3d getCenter() const override { return center; } 

    void buildVoltage() override {
        voltage = -Einc->amplitude
            * conj(getIntegratedPlaneWave(Einc->wavevec).dot(Einc->pol)); // Hermitian dot!
    }

    vec3cd getRadAlongDir(const vec3d& X, const vec3d& kvec) const override {
        return exp(iu*kvec.dot(X)) * getIntegratedPlaneWave(kvec).conjugate();
    }

    vec3cd getFarAlongDir(const vec3d& krhat) const override {
        return getIntegratedPlaneWave(krhat).conjugate();
    }

    void buildSubRWGs();

    vec3cd getIntegratedPlaneWave(const vec3d&, bool = 0) const;

    cmplx getIntegratedRad(const std::shared_ptr<Source>) const override;

private:
    std::array<std::shared_ptr<Triangle>,2> tris;
    std::array<vec3d,2> Xpm;  // non-common vertices
    // std::array<int, 2> idxpm; // global indices of non-common vertices

    int idx0, idx1; // global indices of common vertices
    vec3d X0, X1;   // common vertices

    vec3d center;   // midpoint of common edge
    double leng;    // length of common edge

    std::array<std::shared_ptr<RWG>,14> subrwgs;
    std::unique_ptr<BC> bc;
};