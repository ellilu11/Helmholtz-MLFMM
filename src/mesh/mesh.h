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
    struct TriToRWG;

    using quadPair = std::pair<vec3d, double>;
    using MomentsEFIE = std::tuple<cmplx, vec3cd, vec3cd, cmplx>;
    using MomentsMFIE = std::tuple<cmplx, vec3cd, vec3cd, vec3cd, cmplx>;
    using intRads = std::vector<std::pair<double, vec3d>>;

    // Global data
    std::vector<vec3d> glVerts;   // list of vertices (including fine)
    std::vector<Triangle> glTris; // list of triangles (including fine)
    PairHashMap<size_t> glPairsToIdx; // map from triangle pair to index in glTriPairs
    
    // std::vector<TriToRWG> triToRWGs; // coarse triangle to RWG mappings

    vec3d rootCenter; 
    double rootLeng;

    // Functions
    void buildRootCoords();

    SrcVec importMesh(const std::filesystem::path&);

    void getScattered(const SrcVec&, const std::filesystem::path&, const std::string&, int);

    // void getSurfCurr();

    void printNormals(const std::string&);
}