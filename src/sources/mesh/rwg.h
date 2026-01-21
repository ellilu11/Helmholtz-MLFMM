#pragma once

#include "../source.h"
#include "triangle.h"

class Mesh::RWG : public Source {
public:
    RWG(std::shared_ptr<Excitation::PlaneWave>, size_t, const vec4i&);

    vec3cd getIntegratedPlaneWave(const vec3d&,bool = 0) const;

    cmplx getIntegratedRad(const std::shared_ptr<Source>) const override;

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
            * conj(getIntegratedPlaneWave(Einc->wavevec).dot(Einc->pol)); // Hermitian dot!
    }

    vec3cd getRadAlongDir(const vec3d& X, const vec3d& kvec) const override {
        return exp(iu*kvec.dot(X)) * getIntegratedPlaneWave(kvec).conjugate();
    }

    vec3cd getFarAlongDir(const vec3d& krhat) const override {
        return getIntegratedPlaneWave(krhat).conjugate();
    }

protected:
    vec2i iTris;    // indices of triangles
    vec2i iVertsC;  // indices of common vertices
    vec2i iVertsNC; // indices of non-common vertices

    double leng;    // length of common edge
};