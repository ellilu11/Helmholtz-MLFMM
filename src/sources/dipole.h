#pragma once

#include <random>
#include "source.h"

class Dipole : public Source {
public:
    Dipole() = default;

    Dipole(std::shared_ptr<PlaneWave>, const vec3d&);

    vec3d getCenter() const override { return pos; }

    void buildRHS() override;

    void buildCurrent() override;

    vec3cd getRadAlongDir(const vec3d&, const vec3d&) const override;

    vec3cd getIncAlongDir(const vec3d&, const vec3d&) const override;

    vec3cd getRadAtPoint(const vec3d&) const override;

    cmplx getIntegratedRad(const std::shared_ptr<Source>) const override;

private:
    vec3d pos;  // position
    vec3d phat; // unit pol. density vector
    double pol; // pol. density magnitude 
};