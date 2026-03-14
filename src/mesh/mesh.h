#pragma once

#include <filesystem>
#include "../maps.h"

namespace Mesh {
    // Constants
    constexpr int nCommonThres = 2;
    // constexpr double distThres = 2.0;

    // Types
    class Triangle;
    class RWG;
    struct TriPair;
    struct TriToRWG;

    using MomentsEFIE = std::tuple<cmplx, vec3cd, vec3cd, cmplx>;
    using MomentsMFIE = std::tuple<cmplx, vec3cd, vec3cd, vec3cd, cmplx>;

    // Global data
    std::vector<vec3d> glVerts;   // list of vertices (including fine)
    std::vector<Triangle> glTris; // list of triangles (including fine)
    PairHashMap<TriPair> glTriPairs; // nearfield triangle pairs
    std::vector<TriToRWG> triToRWGs; // coarse triangle to RWG mappings
    size_t nverts;                // number of coarse mesh vertices
    size_t ntris;                 // number of coarse mesh triangles

    vec3d rootCenter; 
    double rootLeng;

    // Functions
    void importVertices(const std::filesystem::path&);

    void importTriangles(const std::filesystem::path&);

    void buildRootCoords();

    SrcVec importRWGs(
        const std::filesystem::path&, std::shared_ptr<Exct::PlaneWave>);

    SrcVec importMesh(
        const std::filesystem::path&, std::shared_ptr<Exct::PlaneWave>);

    void printScattered(const SrcVec&, const std::filesystem::path&, const std::string&, int);

    // void printSurfCurr();

    void printNormals(const std::string&);
}