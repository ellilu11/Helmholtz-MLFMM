#pragma once

#include <filesystem>

namespace Mesh {
    // Types
    class Triangle;
    class RWG;
    struct TriPair;
    struct TriToRWG;

    // Coarse mesh data
    std::vector<vec3d> glVerts;   // list of vertices (including fine)
    std::vector<Triangle> glTris; // list of triangles (including fine)
    PairHashMap<TriPair> glTriPairs; // nearfield triangle pairs
    std::vector<TriToRWG> triToRWGs; // coarse triangle to RWG mappings
    size_t nverts;                // number of coarse mesh vertices
    size_t ntris;                 // number of coarse mesh triangles

    // Functions
    void importVertices(const std::filesystem::path&);

    void importTriangles(const std::filesystem::path&);

    SrcVec importRWGs(
        const std::filesystem::path&,
        const std::shared_ptr<Exc::PlaneWave>);

    SrcVec importMesh(
        const std::filesystem::path&,
        const std::filesystem::path&,
        const std::filesystem::path&,
        std::shared_ptr<Exc::PlaneWave>);

    void printScattered(const SrcVec&, const std::string&, int);

    // void evaluateJ();

    // void printRWGs(const SrcVec& rwgs, const std::string&);
}