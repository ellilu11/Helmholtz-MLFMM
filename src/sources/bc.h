#pragma once

#include "rwg.h"

class BC {
public:
    BC(RWG*);

private:
    // std::array<std::shared_ptr<RWG>, 10> subrwgs;
    std::vector<std::shared_ptr<RWG>> subrwgs;

    // vec2i iVerts; // global indices of common vertices
    vec3d X0;   // 1st common vertex
    vec3d X1;   // 2nd common vertex

    vec3d center; // midpoint of common edge
};