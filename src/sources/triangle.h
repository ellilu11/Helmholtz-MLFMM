#pragma once

#include <filesystem>
#include <ranges>
#include "maps.h"

class Triangle {
public:
    friend class RWG;
    friend class SrcRWG;
    friend class SubRWG;
    friend class BC;

    static void importVertices(const std::filesystem::path& fpath);

    static void importTriangles(
        const std::filesystem::path&,
        const std::filesystem::path&, 
        Precision);

    Triangle(const vec3i&, int);

    static void buildSubtris();

    static void buildSubtriMaps();

    static void buildQuadCoeffs();

    void buildQuads(const std::array<vec3d,3>&);

    std::array<vec3d,3> getVerts() const {
        return { glVerts[iVerts[0]], glVerts[iVerts[1]], glVerts[iVerts[2]] };
    }

    //vec3d getCenter() const {
    //    return (glVerts[iVerts[0]] + glVerts[iVerts[1]] + glVerts[iVerts[2]]) / 3.0;
    //}

    int getBaseIdx() const { return (iTri - NTris) / 6; }

    std::vector<quadPair> getQuads() const { return quads; }

    static int getNumQuads() { return numQuads; }

    // bool isAdjacent(const std::shared_ptr<Triangle>&);

public:
    // TODO: Place in Mesh::
    static std::vector<vec3d> glVerts;     // indices to all vertices
    static std::vector<Triangle> glTris;   // indices to all triangles
    static PairHashMap<int> glEdgeToMid;   // coarse edges to midpoint indices
    static PairHashMap<vec2i> glEdgeToTri; // fine edges to subtri indices 

    static size_t NVerts; // number of coarse mesh vertices
    static size_t NTris;  // number of coarse mesh triangles

private:
    static std::vector<quadPair> quadCoeffs;
    static Precision quadPrec;
    static int numQuads;

    std::vector<quadPair> quads;

    vec3i iVerts; // global indices of vertices
    int iCenter;  // global index of center
    int iTri;     // index in glTris
    // int iBase;    // index of parent coarse triangle (TODO: Move to derived Subtri)

    vec3d center;           // barycentric center
    std::array<vec3d,3> Ds; // edges (Ds[i] = Xs[i+1] - Xs[i])

    vec3d nhat;             // surface normal unit vector
    // double alpha;        // angle between 0th and 2nd edges
};