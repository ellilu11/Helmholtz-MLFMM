#pragma once

#include <ranges>
#include "mesh.h"

class Mesh::Triangle {
public:
    friend class RWG;
    friend class SrcRWG;
    friend class SubRWG;
    friend class BC;

    static void buildQuadCoeffs(Precision);

    Triangle(const vec3i&, int);

    int getNumCommonVerts(const Triangle&) const;

    void buildSelfIntegrated();

    std::pair<double,vec3d> getNearIntegrated(const vec3d&, bool = 0) const;

    std::pair<cmplx,vec3cd> getPlaneWaveIntegrated(const vec3d&) const;

    // cmplx getDuffyIntegrated(const vec3d&, const vec3d&, const vec3d&) const;

    void buildTriQuads();

    static void buildRadMoments();

    static void refineVertices();

    static void refineTriangles();

    static void buildEdgeToTri();

    std::array<vec3d,3> getVerts() const {
        return { glVerts[iVerts[0]], glVerts[iVerts[1]], glVerts[iVerts[2]] };
    }

    std::vector<quadPair<vec3d>> getQuads() const { return triQuads; }

    static int getNumQuads() { return numQuads; }

    vec3d proj(const vec3d& X) const { return X - (nhat.dot(X))*nhat; }

private:
    static int numQuads;
    static std::vector<quadPair<vec3d>> quadCoeffs;
    static std::vector<quadPair<double>> linQuads; // linear quadrature nodes and weights

    std::vector<quadPair<vec3d>> triQuads; // triangle quadrature nodes and weights

    int iTri;     // index in glTris
    vec3i iVerts; // indices of vertices
    int iCenter;  // index of center

    vec3d center;           // barycentric center 
    std::array<vec3d,3> Ds; // edge displacements (Ds[i] = Xs[i+1] - Xs[i])
    vec3d nhat;             // surface normal unit vector
    double area;            // area

    // Integral quantities
    std::array<double,4> selfInts;
};