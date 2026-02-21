#pragma once

#include "../source.h"
#include "triangle.h"

extern int numNearTriPairs;

class Mesh::RWG : public Source {
public:
    RWG(std::shared_ptr<Exc::PlaneWave>, size_t, const vec4i&);

    static void refineRWGs();

    vec3cd getPlaneWaveIntegrated(const vec3d&, bool = 0) const;

    cmplx getIntegratedRad(const std::shared_ptr<Source>) const override;

    std::array<int,2> getTrisIdx() const { return iTris; }

    std::array<Triangle,2> getTris() const {
        return { glTris[iTris[0]], glTris[iTris[1]] };
    }

    std::array<vec3d,2> getVertsC() const {
        return { glVerts[iVertsC[0]], glVerts[iVertsC[1]] };
    }

    std::array<vec3d,2> getVertsNC() const {
        return { glVerts[iVertsNC[0]], glVerts[iVertsNC[1]] };
    }

    vec3d getCenter() const {
        const auto& verts = getVertsC();
        return (verts[0] + verts[1]) / 2.0;
    }

    double getLeng() const { return leng; }

    void buildVoltage() override {
        voltage = -Einc->amplitude
            * conj(getPlaneWaveIntegrated(Einc->wavevec).dot(Einc->pol)); // Hermitian dot!
    }

    vec3cd getRadAlongDir(const vec3d& X, const vec3d& kvec) const override {
        return exp(iu*kvec.dot(X)) * getPlaneWaveIntegrated(kvec).conjugate();
    }

    vec3cd getFarAlongDir(const vec3d& krhat) const override {
        return getPlaneWaveIntegrated(krhat).conjugate();
    }

protected:
    std::array<int,2> iTris;    // indices of triangles
    std::array<int,2> iVertsC;  // indices of common vertices
    std::array<int,2> iVertsNC; // indices of non-common vertices

    double leng; // length of common edge
};