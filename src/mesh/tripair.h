#pragma once

#include "mesh.h"
#include "../types.h"

struct Mesh::TriPair {

public :
    TriPair(pair2i);

    std::pair<Triangle,Triangle> getTriPair() const {
        return { glTris[iTris.first], glTris[iTris.second] };
    }

private :
    void buildNumCommon();

    void buildMomentsEFIE();

    void buildIntegratedInvR();

public : 
    std::tuple<cmplx, vec3cd, vec3cd, cmplx> momentsEFIE;
    std::vector<std::pair<double,vec3d>> integratedInvR; 
    std::vector<std::pair<double, vec3d>> integratedInvR2; // symmetric case

    pair2i iTris; // indices of triangles
    int nCommon;  // number of common vertices
};