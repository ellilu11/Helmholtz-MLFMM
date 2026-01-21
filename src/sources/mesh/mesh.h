#pragma once

#include <filesystem>

namespace Mesh {
    // Mesh classes
    class Triangle;
    class RWG;
    class SrcRWG;
    class SubRWG;
    class BC;

    // Global mesh data
    std::vector<vec3d> glVerts;     // list of vertices (including refined)
    std::vector<Triangle> glTris;   // list of triangles (including refined)
    // SrcVec glRwgs;               // list of rwgs
    std::vector<SubRWG> glSubrwgs;  // list of subrwgs
    size_t nverts;                  // number of coarse mesh vertices
    size_t ntris;                   // number of coarse mesh triangles

    // Mesh refinement maps
    PairHashMap<int> glEdgeToMid;   // coarse edges to midpoint indices
    PairHashMap<vec2i> glEdgeToTri; // fine edges to subtri indices 
    PairHashMap<int> glEdgeToSub;   // edges to subrwg indices
    // indices of subrwg having vert as common vert
    std::vector<std::vector<SubRWG>> vertsToSubrwgs; 

    // Mesh import/refine functions
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