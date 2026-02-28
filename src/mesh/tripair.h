#pragma once

#include "mesh.h"
#include "../types.h"

class Mesh::TriPair {
    friend RWG;

public :
    TriPair(pair2i iTris) 
        : iTris(iTris)
    {
        buildNumCommon();
        buildRadMoments();
        if (nCommon >= 2) buildIntegratedInvR();
    }

    std::pair<Triangle,Triangle> getTriPair() const {
        return { glTris[iTris.first], glTris[iTris.second] };
    }

    double getDoubleIntegratedInvR(const vec3d&, const vec3d&) const;

    int getNumCommonVerts() const { return nCommon; }

private :
    void buildNumCommon();

    void buildRadMoments();

    void buildIntegratedInvR();

private : 
    QuadMoments radMoments; // quadrature moments for non-singular double integrals
    std::vector<std::pair<double,vec3d>> integratedInvR; // integrated 1/R term

    pair2i iTris; // indices of triangles
    int nCommon;  // number of common vertices
};