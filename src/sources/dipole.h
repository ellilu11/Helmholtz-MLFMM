#pragma once

#include <random>
#include "source.h"

class Dipole : public Source {
public:
    Dipole() = default;

    Dipole(const std::shared_ptr<PlaneWave>&, const vec3d& X) // TODO: Pass in pol. density vector
        : Source(Einc), pos(X), phat(vec3d(1,0,0)), pol(1.0) { 
    };

    vec3d getCenter() const override { return pos; }

    void buildRHS() override;

    void buildCurrent() override;

    vec3cd getRadAlongDir(const vec3d&, const vec3d&) const override;

    vec3cd getIncAlongDir(const vec3d&, const vec3d&) const override;

    vec3cd getRad(const vec3d&, double) const override;

    cmplx getIntegratedRad(const std::shared_ptr<Source>, double) const override;

private:
    vec3d pos;  // position
    vec3d phat; // unit pol. density vector
    double pol; // pol. density magnitude 
};