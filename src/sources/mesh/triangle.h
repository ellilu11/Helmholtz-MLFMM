#pragma once

#include <ranges>
#include "mesh.h"

class Mesh::Triangle {
public:
    friend class RWG;
    friend class SrcRWG;
    friend class SubRWG;
    friend class BC;

    Triangle(const vec3i&, int);

    static void refineVertices();

    static void refineTriangles();

    static void buildEdgeToTri();

    static void buildQuadCoeffs(Precision);

    void buildQuads(const std::array<vec3d,3>&);

    std::array<vec3d,3> getVerts() const {
        return { glVerts[iVerts[0]], glVerts[iVerts[1]], glVerts[iVerts[2]] };
    }

    //vec3d getCenter() const {
    //    return (glVerts[iVerts[0]] + glVerts[iVerts[1]] + glVerts[iVerts[2]]) / 3.0;
    //}

    std::vector<quadPair> getQuads() const { return quads; }

    static int getNumQuads() { return numQuads; }

private:
    static std::vector<quadPair> quadCoeffs;
    static Precision quadPrec;
    static int numQuads;

    std::vector<quadPair> quads;

    int iTri;     // index in glTris
    vec3i iVerts; // indices of vertices
    int iCenter;  // index of center

    // Store or compute these on the fly?
    vec3d center;           // barycentric center 
    std::array<vec3d,3> Ds; // edge displacements (Ds[i] = Xs[i+1] - Xs[i])
    vec3d nhat;             // surface normal unit vector
    double area;            // area of triangle
};