#pragma once

#include <ranges>
#include "mesh.h"

class Mesh::Triangle {
    friend class TriPair;
    friend class TriMoments;
    friend class RWG;

public:
    static void buildQuadCoeffs(Precision);

    Triangle(const vec3i&, int);

    void reverseOrient();

    void buildSelfIntegratedInvR();

    double getDoubleSelfIntegratedInvR(const vec3d&, const vec3d&) const;

    std::pair<cmplx, vec3cd> getIntegratedPlaneWave(const vec3d&) const;

    std::pair<double,vec3d> getIntegratedInvR(const vec3d&, bool = false) const;

    std::pair<double, vec3d> getIntegratedInvRcubed(const vec3d&, bool = false) const;

    double getSingularEFIE(
        const Triangle&, const TriPair&, const vec3d&, const vec3d&) const;

    double getSingularMFIE(
        const Triangle&, const TriPair&, const vec3d&, const vec3d&) const;

    // cmplx getSurfaceCurrent() const;

    void buildTriQuads();

    static void refineVertices();

    static void refineTriangles();

    static void buildEdgeToTri();

    vec3d proj(const vec3d& X) const { return X - (nhat.dot(X))*nhat; }

    std::array<vec3d,3> getVerts() const {
        return { glVerts[iVerts[0]], glVerts[iVerts[1]], glVerts[iVerts[2]] };
    }

    std::vector<quadPair> getQuads() const { return quads; }

    int getIdx() const { return iTri; }

    vec3d getCenter() const { return center; }

    vec3d getNormal() const { return nhat; }

private:
    static std::vector<quadPair> quadCoeffs;

    // Integral quantities
    std::vector<quadPair> quads; // triangle quadrature nodes and weights
    std::array<double, 4> selfInts;

    vec3d center;   // barycentric center 
    vec3d nhat;     // surface normal unit vector
    double area;    

    vec3i iVerts;   // indices of vertices
    // int iCenter; // index of center
    const int iTri; // index in glTris
};

struct Mesh::TriToRWG {
    std::vector<int> iRWGs; // indices of RWGs containing this triangle
    std::vector<int> isMinus; // whether this triangle is positive or negative in each RWG
};