#pragma once

class Triangle;

using TriVec = std::vector<std::shared_ptr<Triangle>>;

class Triangle {
public:
    friend class RWG;

    Triangle()
        : vIdx(vec3i::Zero()),
          Xs( {zeroVec,zeroVec,zeroVec})
          {};

    Triangle(
        const vec3i&,
        const std::vector<vec3d>&,
        const Precision);

    std::vector<quadPair> getQuads() { return quads; }

    static int getNumQuads(Precision);

    static void buildQuadCoeffs(Precision);

    void buildQuads(Precision);

    // bool isAdjacent(const std::shared_ptr<Triangle>&);

private:
    static std::vector<quadPair> quadCoeffs;
    std::vector<quadPair> quads;

    vec3i vIdx;
    std::array<vec3d,3> Xs; // vertices
    std::array<vec3d,3> Ds; // Ds[i] = Xs[i+1] - Xs[i]
    vec3d nhat; // surface normal unit vector
};