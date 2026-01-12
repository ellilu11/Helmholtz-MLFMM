#pragma once

class Triangle;

using TriVec = std::vector<std::shared_ptr<Triangle>>;

class Triangle {
public:
    friend class RWG;

    //Triangle()
    //    : vIdx(vec3i::Zero()),
    //      Xs( {zeroVec,zeroVec,zeroVec})
    //      {};

    Triangle(const vec3i&, const std::vector<vec3d>&, Precision);

    Triangle(int, const vec3d&, const vec3d&, const vec3d&, Precision);

    void buildTriangle(Precision);

    std::vector<quadPair> getQuads() { return quads; }

    static void buildQuadCoeffs(Precision);

    static int getNumQuads() { return numQuads; }

    void buildQuads(Precision);

    void findRefinedTris(Precision);

    // bool isAdjacent(const std::shared_ptr<Triangle>&);

private:
    static std::vector<quadPair> quadCoeffs;
    static int numQuads;

    static TriVec refinedTris;

    std::vector<quadPair> quads;

    vec3i vIdx;             // global indices of vertices
    std::array<vec3d,3> Xs; // vertices
    std::array<vec3d,3> Ds; // Ds[i] = Xs[i+1] - Xs[i]
    vec3d center;           // barycentric center
    vec3d nhat;             // surface normal unit vector
};