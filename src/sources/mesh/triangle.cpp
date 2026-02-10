#include "triangle.h"

int Mesh::Triangle::numQuads;
std::vector<quadPair<vec3d>> Mesh::Triangle::quadCoeffs;
std::vector<quadPair<double>> Mesh::Triangle::linQuads;

/* Triangle(iVerts)
 * Construct triangle from global indices of vertices
 * iVerts : global indices of vertices
 */
Mesh::Triangle::Triangle(const vec3i& iVerts, int iTri)
    : iVerts(iVerts), iTri(iTri)
{
    assert(iVerts.maxCoeff() < glVerts.size());

    const auto Xs = getVerts();
    center = (Xs[0] + Xs[1] + Xs[2]) / 3.0;

    for (int i = 0; i < 3; ++i) Ds[i] = Xs[(i+1)%3] - Xs[i];
    nhat = (Ds[0].cross(Ds[1])).normalized();

    // If nhat is pointing inward, reverse it
    // Assume a closed, star-shaped mesh centered at and enclosing the origin
    if (center.dot(nhat) < 0.0) nhat *= -1.0;

    area = (Ds[0].cross(-Ds[2])).norm() / 2.0;

    buildTriQuads(Xs);

    // std::cout << "Built triangle #" << iTri << " with area " << area << '\n';
}

// Return global indices of vertices shared by this and other triangle
intVec Mesh::Triangle::getCommonVerts(const Triangle& other) const {
    // if (iTri == other.iTri) return ;

    intVec iCommons;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (iVerts[i] == other.iVerts[j])
                iCommons.push_back(iVerts[i]);

    assert(iCommons.size() <= 3);
    return iCommons;
}

cmplx Mesh::Triangle::getAdjacentIntegrated(
    const vec3d& src, const vec3d& obsNC, const vec3d& srcNC) const
{
    cmplx sclRad = 0.0, vecRad = 0.0;

    for (int i = 0; i < 3; ++i) {
        const auto& V0 = glVerts[iVerts[i]];
        const auto& V1 = glVerts[iVerts[(i+1)%3]];

    }
    

    return sclRad + vecRad;
}

cmplx Mesh::Triangle::getDuffyIntegrated(
    const vec3d& src, const vec3d& obsNC, const vec3d& srcNC) const 
{
    const double k = FMM::wavenum;

    cmplx rad = 0.0;
    // double areaSum = 0.0;
    for (int i = 0; i < 3; ++i) {
        const auto& V0 = glVerts[iVerts[i]];
        const auto& V1 = glVerts[iVerts[(i+1)%3]];
        const double area = (V0-src).cross(V1-src).norm()/2.0; // cache?
        if (Math::approxZero(area)) continue;

        // areaSum += area;

        cmplx triRad = 0.0;
        for (const auto& [u, uweight] : linQuads) {
            for (const auto& [v, vweight] : linQuads) {
                const vec3d& obs = src + u*(V0 - src) + u*v*(V1 - V0);

                const double dr = (obs-src).norm();
                const double alpha = dr / u;

                triRad += 
                    ((obs-obsNC).dot(src-srcNC) - 4.0 / (k*k))
                    * exp(iu*k*dr) / alpha // u * Math::helmholtzG(dr, k)
                    * uweight * vweight;
            }
        }

        rad += area * triRad;
    }

    // assert(Math::approxZero(areaSum-area));
    // std::cout << areaSum << ' ' << area << '\n';

    return rad;
}

void Mesh::Triangle::buildRadMoments() {
    /* TODO: Only loop over triangles in nearfield of each other
    for (size_t iObs = 0; iObs < glTris.size(); ++iObs) {
        for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) {
            const auto& obsTri = glTris[iObs];
            const auto& srcTri = glTris[iSrc];

            const auto key = makeUnordered(obsTri.iTri, srcTri.iTri);

            QuadMoments moments = { 0.0, vec3cd::Zero(), vec3cd::Zero(), 0.0 };

            for (const auto& [obs, obsWeight] : obsTri.quads) {
                for (const auto& [src, srcWeight] : srcTri.quads) {
                    const cmplx coeff = 
                        Math::helmholtzG((obs-src).norm(), FMM::wavenum) * obsWeight * srcWeight;
                    get<0>(moments) += obs.dot(src) * coeff;
                    get<1>(moments) += src * coeff;
                    get<2>(moments) += obs * coeff;
                    get<3>(moments) += coeff;
                }
            }

            glRadMoments.emplace(key, moments);
        }
    }
    */
}

// Refine vertices: add centers and midpoints of coarse tris
void Mesh::Triangle::refineVertices() {
    int iVerts = nverts;

    for (auto& tri : glTris) {
        // Add center
        tri.iCenter = iVerts++; // glVerts.size();
        glVerts.push_back(tri.center);

        // Add midpoints of edges
        for (int i = 0; i < 3; ++i) {
            const int idx0 = tri.iVerts[i], idx1 = tri.iVerts[(i+1)%3];
            const vec3d mid = (glVerts[idx0] + glVerts[idx1]) / 2.0;
            const auto& edge = makeUnordered(idx0, idx1);

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
                const auto& iVerts = (!j ?
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
            const auto& edge = makeUnordered(idx0, idx1);

            if (fineEdgeToTri.find(edge) == fineEdgeToTri.end())
                fineEdgeToTri.emplace(edge, vec2i(tri.iTri, -1));
            else
                fineEdgeToTri[edge][1] = tri.iTri;
        }
    }
}

