#pragma once

#include <filesystem>
#include "maps.h"

class Triangle;

using TriVec = std::vector<std::shared_ptr<Triangle>>;
using TriArr6 = std::array<std::shared_ptr<Triangle>,6>;

class Triangle {
public:
    friend class RWG;
    friend class SrcRWG;
    friend class SubRWG;
    friend class BC;

    static void importVertices(const std::filesystem::path& fpath);

    static TriVec importTriangles(
        const std::filesystem::path&,
        const std::filesystem::path&, 
        Precision);

    Triangle(const vec3i&);

    Triangle(int, int, int);

    static void refineVertices(const TriVec&);

    TriArr6 getSubtris();

    static void buildQuadCoeffs();

    void buildQuads(const std::array<vec3d,3>&);

    std::array<vec3d,3> getVerts() const {
        return { glVerts[glIdxs[0]], glVerts[glIdxs[1]], glVerts[glIdxs[2]] };
    }

    //vec3d getCenter() const {
    //    return (glVerts[glIdxs[0]] + glVerts[glIdxs[1]] + glVerts[glIdxs[2]]) / 3.0;
    //}

    std::vector<quadPair> getQuads() { return quads; }

    static int getNumQuads() { return numQuads; }

    // bool isAdjacent(const std::shared_ptr<Triangle>&);

public:
    static std::vector<vec3d> glVerts;   // global list of fine vertices
    static PairHashMap<int> glEdgeToIdx; // global list of indices of edges

private:
    static std::vector<quadPair> quadCoeffs;
    static Precision quadPrec;
    static int numQuads;

    std::vector<quadPair> quads;

    vec3i glIdxs;           // global indices of vertices
    int glIdxCenter;        // global index of barycentric center

    // std::array<vec3d,3> Xs; // vertices
    vec3d center;           // barycentric center
    std::array<vec3d,3> Ds; // edges (Ds[i] = Xs[i+1] - Xs[i])

    vec3d nhat;             // surface normal unit vector
    // double alpha;        // angle between 0th and 2nd edges
};