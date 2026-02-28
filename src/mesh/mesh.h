#pragma once

#include <filesystem>

namespace Mesh {
    // Constants
    constexpr std::array<double,5> rcoeffs =
        { -1.0/6.0, 1.0/6.0, -1.0/3.0, 1.0/3.0, 1.0 };

    // Types
    class Triangle;
    class RWG;
    class SrcRWG;
    class SubRWG;
    class BC;
    class TriPair;
    struct TriToRWG;

    using QuadMoments = std::tuple<cmplx, vec3cd, vec3cd, cmplx>;

    // Coarse mesh data
    std::vector<vec3d> glVerts;   // list of vertices (including fine)
    std::vector<Triangle> glTris; // list of triangles (including fine)
    PairHashMap<TriPair> glTriPairs; // nearfield triangle pairs
    std::vector<TriToRWG> triToRWGs; // coarse triangle to RWG mappings
    size_t nverts;                // number of coarse mesh vertices
    size_t ntris;                 // number of coarse mesh triangles

    // Fine mesh data
    std::vector<SubRWG> glSubrwgs;  // list of fine RWGs
    std::vector<double> massCoeffs; // mass coefficients between fine RWGs
    PairHashMap<int> idxMassCoeffs; // subRWG pairs to mass coeff indices

    // Fine mesh maps
    PairHashMap<int> edgeToMid;       // coarse edges to midpoint indices
    PairHashMap<vec2i> fineEdgeToTri; // fine edges to subtri indices 
    PairHashMap<int> fineEdgeToSub;   // fine edges to subrwg indices
    std::vector<std::vector<int>> triToSubs;    // indices of subrwgs containing tri
    std::vector<std::vector<int>> vertToSubs;   // indices of subrwgs containing vert

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

    void refineMesh(const SrcVec&);

    void printScattered(const SrcVec&, const std::string&, int);

    // void evaluateJ();

    // void printRWGs(const SrcVec& rwgs, const std::string&);
}