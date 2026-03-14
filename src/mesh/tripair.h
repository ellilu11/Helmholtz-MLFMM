#pragma once

#include "triangle.h"
#include "../types/types.h"

struct Mesh::TriPair {

public:
    TriPair(pair2i);

    std::pair<Triangle, Triangle> getTriPair() const {
        return { glTris[iTris.first], glTris[iTris.second] };
    }

private:
    void buildNumCommon();

    void buildMomentsEFIE();

    void buildMomentsMFIE();

    void buildIntegratedInvR();

    void buildIntegratedInvRcubed();

public:
    std::unique_ptr<MomentsEFIE> momentsEFIE;
    std::unique_ptr<MomentsMFIE> momentsMFIE;
    std::unique_ptr<MomentsMFIE> momentsMFIE2; // symmetric case

    std::vector<std::pair<double, vec3d>> integratedInvR;
    std::vector<std::pair<double, vec3d>> integratedInvR2; // symmetric case

    std::vector<std::pair<double, vec3d>> integratedInvRcubed;
    std::vector<std::pair<double, vec3d>> integratedInvRcubed2; // symmetric case

    pair2i iTris; // indices of triangles
    double dist;  // distance between triangle centers
    int nCommon;  // number of common vertices
};