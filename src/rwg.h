#pragma once

#include <algorithm>
#include <filesystem>
#include <random>
#include "math.h"
#include "source.h"
#include "triangle.h"

class RWG;

using RWGVec = std::vector<std::shared_ptr<RWG>>;

class RWG {
public:
    RWG(const Eigen::Vector4i&, 
        const std::vector<vec3d>&,
        const TriVec&,
        const std::shared_ptr<Source>);

    std::shared_ptr<Triangle> getTriPlus() const { return triPlus; }

    std::shared_ptr<Triangle> getTriMinus() const { return triMinus; }

    vec3d getCenter() const { return vCenter; }

    vec3d getVplus() const { return vPlus; }

    vec3d getVminus() const { return vMinus; }

    double getLeng() const { return leng; }

    double getCurrent() const { return current; }

    cmplx getSol() const { return sol; }

    void addToSol(cmplx sol_) { sol += sol_; }

    void resetSol() { sol = 0.0; }

    void buildRHS();

    void buildCurrent();

    vec3cd getRadAlongDir(const vec3d&, const vec3d&) const;

    vec3cd getIncAlongDir(const vec3d&, const vec3d&) const;

    vec3cd getRad(const vec3d&, double) const;

    cmplx getIntegratedRad(const std::shared_ptr<RWG>, double) const;

private:
    vec3d v0;
    vec3d v1;
    std::shared_ptr<Triangle> triPlus;
    std::shared_ptr<Triangle> triMinus;

    vec3d vCenter;
    vec3d vPlus;
    vec3d vMinus;

    vec3d nhat; // unit vector normal to RWG surface

    double leng;

    std::shared_ptr<Source> Einc;
    cmplx rhs;
    double current;
    cmplx sol;

};