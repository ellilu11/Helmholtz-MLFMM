#pragma once

#include "../exct.h"
#include "../maps.h"

namespace Mesh {
    // Constants
    constexpr int nCommonThres = 2;
    // constexpr double distThres = 2.0;

    // Types
    class Triangle;
    class RWG;
    struct TriPair;
    struct TriMoments;
    struct TriToRWG;

    using quadPair = std::pair<vec3d, double>;
    using MomentsEFIE = std::tuple<cmplx, vec3cd, vec3cd, cmplx>;
    using MomentsMFIE = std::tuple<cmplx, vec3cd, vec3cd, vec3cd, cmplx>;

    // Global data
    std::vector<vec3d> glVerts;   // list of vertices (including fine)
    std::vector<Triangle> glTris; // list of triangles (including fine)
    PairHashMap<TriPair> glTriPairs; // nearfield triangle pairs
    std::vector<MomentsEFIE> glMomentsEFIE;
    std::vector<MomentsMFIE> glMomentsMFIE;
    std::vector<MomentsMFIE> glMomentsMFIE2; // symmetric case

    std::vector<TriToRWG> triToRWGs; // coarse triangle to RWG mappings

    vec3d rootCenter; 
    double rootLeng;

    // Functions
    void buildRootCoords();

    SrcVec importMesh(const std::filesystem::path&);

    void getScattered(const SrcVec&, const std::filesystem::path&, const std::string&, int);

    // void getSurfCurr();

    void printNormals(const std::string&);
}