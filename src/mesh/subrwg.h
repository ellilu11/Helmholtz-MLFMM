#pragma once

#include "rwg.h"

class Mesh::SubRWG final : public RWG {
    friend class BC;

public:
    SubRWG(int, const vec4i&);

    static void buildVertsToSubRWGs(int);

    static void buildMassCoeffs();

    double getMassCoeff(const SubRWG&) const;

    void setOriented(int, const vec3d&, const vec3d&);

private:
    realVec coeffs;

    // global index of vertex (if it exists) in coarse mesh contributing to BC
    std::optional<int> iVertsCoarse; 

    double oriented;
};