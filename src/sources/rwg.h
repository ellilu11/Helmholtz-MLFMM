#pragma once

#include "source.h"
#include "triangle.h"

class RWG : public Source {
public:
    RWG(const std::shared_ptr<PlaneWave>&,
        const Eigen::Vector4i&, 
        const std::vector<vec3d>&,
        const TriVec&);

    std::shared_ptr<Triangle> getTriPlus() const { return triPlus; }

    std::shared_ptr<Triangle> getTriMinus() const { return triMinus; }

    vec3d getVplus() const { return vPlus; }

    vec3d getVminus() const { return vMinus; }

    double getLeng() const { return leng; }

    void buildRHS() override;

    void buildCurrent() override;

    vec3d getCenter() const override { return center; } 

    vec3cd getRadAlongDir(const vec3d&, const vec3d&) const override;

    vec3cd getIncAlongDir(const vec3d&, const vec3d&) const override;

    vec3cd getRad(const vec3d&, double) const override;

    cmplx getIntegratedRad(const std::shared_ptr<Source>, double) const override;

private:
    vec3d v0;
    vec3d v1;
    vec3d center;

    std::shared_ptr<Triangle> triPlus;
    std::shared_ptr<Triangle> triMinus;
    vec3d vPlus;
    vec3d vMinus;

    vec3d nhat; // unit vector normal to RWG surface

    double leng;
};