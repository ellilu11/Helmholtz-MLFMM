#include "../triangle.h"

/* If nhat is pointing inward, reverse it
// Assume a closed, star-shaped mesh centered at and enclosing the origin
void Mesh::Triangle::fixNhat() {}
    if (center.dot(nhat) < 0.0) {
        nhat *= -1.0;
        std::swap(this->iVerts[0], this->iVerts[2]); // Swap verts per RHR orientation
}
*/

// Refine vertices: add centers and midpoints of coarse tris
void Mesh::Triangle::refineVertices() {
    int iVerts = nverts;

    for (auto& tri : glTris) {
        // Add center
        tri.iCenter = iVerts++; // glVerts.size();
        glVerts.push_back(tri.center);

        // Add midpoints of edges
        for (int i = 0; i < 3; ++i) {
            int idx0 = tri.iVerts[i], idx1 = tri.iVerts[(i+1)%3];
            vec3d mid = (glVerts[idx0] + glVerts[idx1]) / 2.0;
            pair2i edge = makeUnordered(idx0, idx1);

            // Check if midpoint already exists
            if (edgeToMid.find(edge) != edgeToMid.end())
                continue;

            edgeToMid.emplace(edge, iVerts++); // glVerts.size();
            glVerts.push_back(mid);
        }
    }

    assert(iVerts == glVerts.size());
}

// Construct all 6 fine tris of all coarse tris and add to global list
void Mesh::Triangle::refineTriangles() {
    std::vector<Triangle> glSubtris;
    glSubtris.reserve(6 * glTris.size());

    int iSubtri = glTris.size();
    for (const auto& tri : glTris) {
        for (int i = 0; i < 3; ++i) {
            const int idx0 = tri.iVerts[i], idx1 = tri.iVerts[(i+1)%3], idxCenter = tri.iCenter;
            const int idxMid = edgeToMid[makeUnordered(idx0, idx1)];

            for (int j = 0; j < 2; ++j) {
                // Assign vertex on coarse mesh as first vertex of subtri
                vec3i iVerts = (!j ?
                    vec3i(idx0, idxMid, idxCenter) :
                    vec3i(idx1, idxCenter, idxMid)
                    );

                glSubtris.emplace_back(iVerts, iSubtri);
                assert(Math::vecEquals(tri.nhat, glSubtris[iSubtri-glTris.size()].nhat)); // Check orientation
                ++iSubtri;
            }
        }
    }

    // Append list of subtris to global list of tris
    glTris.insert(glTris.end(), glSubtris.begin(), glSubtris.end());
}

// TODO: Move into refineTriangles()
// Build edge to subtri map
void Mesh::Triangle::buildEdgeToTri() {
    const auto glSubtris = glTris | std::views::drop(ntris);

    for (const auto& tri : glSubtris) {

        for (int i = 0; i < 3; ++i) {
            const int idx0 = tri.iVerts[i], idx1 = tri.iVerts[(i+1)%3];
            pair2i edge = makeUnordered(idx0, idx1);

            if (fineEdgeToTri.find(edge) == fineEdgeToTri.end())
                fineEdgeToTri.emplace(edge, vec2i(tri.iTri, -1));
            else
                fineEdgeToTri[edge][1] = tri.iTri;
        }
    }
}