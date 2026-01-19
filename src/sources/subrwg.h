#pragma once

#include "rwg.h"

class SubRWG final : public RWG {
    friend class BC;

public:
    static void buildSubRWGs();

    SubRWG(const vec4i&, int);

    static void buildVertsToSubrwgs(int);

    void setOriented(const vec3d&, const vec3d&, const vec3d&);

private:
    pair2i edge;  // global indices of common edge vertices

    // global index of vertex (if it exists) in coarse mesh contributing to BC
    std::optional<int> iVertsCoarse; 

    double oriented;
};