#pragma once

#include "../source.h"
#include "triangle.h"
#include "tripair.h"

class Mesh::RWG : public Source {
public:
    RWG(std::shared_ptr<Exct::PlaneWave>, size_t, const vec4i&);

    void buildVoltage() override;

    std::pair<vec3cd,vec3cd> getIntegratedPlaneWave(const vec3d&, bool = 0) const;

    cmplx getIntegratedEFIE(const std::shared_ptr<Source>) const override;

    cmplx getIntegratedMFIE(const std::shared_ptr<Source>) const override;

    double getIntegratedMass(const std::shared_ptr<Source>) const override;

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

    std::array<std::pair<Triangle, vec3d>, 2> getTrisAndVerts() const {
        return { std::make_pair(glTris[iTris[0]], glVerts[iVertsNC[0]]),
                 std::make_pair(glTris[iTris[1]], glVerts[iVertsNC[1]]) };
    }

    vec3d getCenter() const {
        auto [v0, v1] = getVertsC();
        return (v0 + v1) / 2.0;
    }

    double getLeng() const { return leng; }

    std::pair<vec3cd,vec3cd> getRadsAlongDir(const vec3d& X, const vec3d& kvec) const override {
        cmplx exp = std::exp(iu*kvec.dot(X)) / (4.0*PI); // apply factor of 1/(4pi)
        auto [rad, radNormal] = getIntegratedPlaneWave(kvec);

        return { exp * rad.conjugate(), exp * radNormal.conjugate() };
    }

    vec3cd getFarAlongDir(const vec3d& krhat) const override {
        return getIntegratedPlaneWave(krhat).first.conjugate() / (4.0*PI); // apply factor of 1/(4pi)
    }

protected:
    std::array<int,2> iTris;    // indices of triangles
    std::array<int,2> iVertsC;  // indices of common vertices
    std::array<int,2> iVertsNC; // indices of non-common vertices

    double leng; // length of common edge
};