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

    Triangle(int idx, const vec3d&, const vec3d&);

    std::vector<quadPair> getQuads() { return quads; }

    static void buildQuadCoeffs();

    static int getNumQuads() { return numQuads; }

    void buildQuads();

    TriArr6 getSubtris(const vec3i&);

    // bool isAdjacent(const std::shared_ptr<Triangle>&);

private:
    static std::vector<vec3d> glVerts;

    static std::vector<quadPair> quadCoeffs;
    static Precision quadPrec;
    static int numQuads;

    std::vector<quadPair> quads;

    vec3i glIdxs;           // global indices of vertices
    // vec3i llIdxs;           // local indices of vertices
    std::array<vec3d,3> Xs; // vertices
    std::array<vec3d,3> Ds; // edges (Ds[i] = Xs[i+1] - Xs[i])
    vec3d center;           // barycentric center
    vec3d nhat;             // surface normal unit vector
};