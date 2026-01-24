#pragma once

#include "rwg.h"

class Mesh::SubRWG final : public RWG {
    friend class BC;

public:
    SubRWG(int, const vec4i&);

    static void buildVertsToSubRWGs(int);

    void setOriented(int, const vec3d&, const vec3d&);

    //std::vector<int> getBases() const { return iBases; }

    //void addBase(int iBase) { iBases.push_back(iBase); }

private:
    // std::vector<int> iBases; // indices of parent RWGs
    realVec coeffs;

    // global index of vertex (if it exists) in coarse mesh contributing to BC
    std::optional<int> iVertsCoarse; 

    double oriented;
};