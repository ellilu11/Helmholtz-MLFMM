#pragma once

#include "source.h"
#include "triangle.h"

class BC final : public Source {
public:
    BC(std::shared_ptr<Excitation::PlaneWave>,
       size_t,
       vec3d, vec3d);

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

    vec3cd getIntegratedPlaneWave(const vec3d&, bool = 0) const;

    cmplx getIntegratedRad(const std::shared_ptr<Source>) const override;

private:
    TriVec tris;

    vec3d X0; // 1st common vertex
    vec3d X1; // 2nd common vertex

    vec3d center; // midpoint of common edge
};