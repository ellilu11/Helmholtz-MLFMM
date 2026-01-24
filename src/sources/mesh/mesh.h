#pragma once

#include <filesystem>

namespace Mesh {
    // Types
    class Triangle;
    class RWG;
    class SrcRWG;
    class SubRWG;
    class BC;

    // Global data
    std::vector<vec3d> glVerts;     // list of vertices (including fine)
    std::vector<Triangle> glTris;   // list of triangles (including fine)
    // SrcVec glRwgs;               // list of rwgs (store globally?)
    std::vector<SubRWG> glSubrwgs;  // list of fine rwgs
    size_t nverts;                  // number of coarse mesh vertices
    size_t ntris;                   // number of coarse mesh triangles

    // Refinement maps
    PairHashMap<int> edgeToMid;       // coarse edges to midpoint indices
    PairHashMap<vec2i> fineEdgeToTri; // fine edges to subtri indices 
    PairHashMap<int> fineEdgeToSub;   // fine edges to subrwg indices
    // std::vector<intVec> vertToTris; // indices of tris containing vert
    std::vector<intVec> vertToSubs; // indices of subrwgs containing vert

    // Functions
    void importVertices(const std::filesystem::path&);

    void importTriangles(const std::filesystem::path&);

    SrcVec importRWGs(
        const std::filesystem::path&,
        const std::shared_ptr<Excitation::PlaneWave>);

    SrcVec importMesh(
        const std::filesystem::path&,
        const std::filesystem::path&,
        const std::filesystem::path&,
        const std::shared_ptr<Excitation::PlaneWave>);

    void refineMesh(const SrcVec&);
}