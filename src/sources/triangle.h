#pragma once

#include <filesystem>

class Triangle;

using TriVec = std::vector<std::shared_ptr<Triangle>>;
using TriArr6 = std::array<std::shared_ptr<Triangle>,6>;

class Triangle {
public:
    friend class RWG;
    friend class BC;

    static void importVertices(const std::filesystem::path& fpath);

    static TriVec importTriangles(
        const std::filesystem::path&,
        const std::filesystem::path&, 
        Precision);

    Triangle(const vec3i&);

    Triangle(const std::array<vec3d,3>&);

    void buildTriangle();

    std::vector<quadPair> getQuads() { return quads; }

    static void buildQuadCoeffs();

    static void buildVertsToSubtris(int);

    static int getNumQuads() { return numQuads; }

    void buildQuads();

    TriArr6 getSubtris(const std::array<vec3d,3>&);

    // bool isAdjacent(const std::shared_ptr<Triangle>&);

private:
    static std::vector<vec3d> glVerts;

    static std::vector<quadPair> quadCoeffs;
    static Precision quadPrec;
    static int numQuads;

    //static TriVec glSubtris; // list of all subtris in mesh
    //static std::vector<std::vector<int>> vertsToSubtris; // list of indices of subtri containing vert

    std::vector<quadPair> quads;

    vec3i iVerts;           // global indices of vertices
    std::array<vec3d,3> Xs; // vertices
    std::array<vec3d,3> Ds; // edges (Ds[i] = Xs[i+1] - Xs[i])
    vec3d center;           // barycentric center
    vec3d nhat;             // surface normal unit vector
};