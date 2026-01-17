#pragma once

#include "rwg.h"

class SubRWG final : public RWG {
    friend class BC;

public:
    SubRWG(std::shared_ptr<Triangle>, std::shared_ptr<Triangle>);

    static void buildVertsToSubrwgs(int);

    void setOriented(const vec3d&, const vec3d&, const vec3d&);

private:
    // global index of vertex (if it exists) in coarse mesh contributing to BC
    std::optional<int> glIdx; 

    double oriented;
};