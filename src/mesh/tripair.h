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

    void buildIntegratedInvR();

    void buildIntegratedInvRcubed();

public:
    std::vector<std::pair<double, vec3d>> intInvR;
    std::vector<std::pair<double, vec3d>> intInvR2; // symmetric case

    std::vector<std::pair<double, vec3d>> intInvRcubed;
    std::vector<std::pair<double, vec3d>> intInvRcubed2; // symmetric case

    pair2i iTris; // indices of triangles
    double dist;  // distance between triangle centers
    int nCommon;  // number of common vertices
    size_t iPair;    // index in glTriPairs
};

struct Mesh::TriMoments {

public:
    TriMoments(size_t);

    void clear() {
        momentsEFIE.clear();
        momentsMFIE.clear();
        momentsMFIE2.clear();
    }

private:
    void buildMomentsEFIE();

    void buildMomentsMFIE();

public:


private:
    size_t nPair; // number of triangle pairs
};