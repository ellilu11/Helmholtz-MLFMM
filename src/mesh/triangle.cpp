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

    const auto& Xs = getVerts();
    center = (Xs[0] + Xs[1] + Xs[2]) / 3.0;

    for (int i = 0; i < 3; ++i) Ds[i] = Xs[(i+1)%3] - Xs[i];
    nhat = (Ds[0].cross(Ds[1])).normalized();

    /* If nhat is pointing inward, reverse it
    // Assume a closed, star-shaped mesh centered at and enclosing the origin
    if (center.dot(nhat) < 0.0) {
        nhat *= -1.0;
        std::swap(this->iVerts[0], this->iVerts[2]); // Swap verts per RHR orientation
    }
    */

    area = (Ds[0].cross(-Ds[2])).norm() / 2.0;

    buildTriQuads();
    buildSelfIntegrated();
    // std::cout << "Built triangle #" << iTri << " with nhat " << nhat << '\n';
}

// Return global indices of vertices shared by this and other triangle
int Mesh::Triangle::getNumCommonVerts(const Triangle& other) const {
    if (iTri == other.iTri) return 3;

    int numVerts = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (iVerts[i] == other.iVerts[j])
                ++numVerts;

    assert(numVerts < 3);
    return numVerts;
}

std::tuple<double,vec3d,vec3d> 
    Mesh::Triangle::getNearIntegrated(const vec3d& obs, bool doNumeric) const
{
    using namespace Math;

    double scaRad = 0.0;
    vec3d vecRad = vec3d::Zero();
    const vec3d& R = projectToPlane(obs);

    if (doNumeric) {
        for (const auto& [node, weight] : triQuads) {
            const double dr = (obs-node).norm();
            scaRad += weight / dr;
            vecRad -= weight * (R-projectToPlane(node)) / dr;
        }

        return std::make_tuple(scaRad, vecRad, R);
    }

    const auto& verts = getVerts();
    double d = std::fabs(nhat.dot(obs - verts[0])), dsq = d*d;

    for (int i = 0; i < 3; ++i) {
        const vec3d& V0 = verts[i];
        const vec3d& V1 = verts[(i+1)%3];

        const vec3d& P0 = projectToPlane(V0) - R;
        const vec3d& P1 = projectToPlane(V1) - R;

        const vec3d& lhat = (P1-P0).normalized();
        const vec3d& uhat = lhat.cross(nhat); 
        double l0 = P0.dot(lhat), l1 = P1.dot(lhat);

        const vec3d& P = P0 - l0*lhat;
        double p = P.norm(), p0 = P0.norm(), p1 = P1.norm();
        const vec3d& phat = (approxZero(p) ? zeroVec : P.normalized()); 

        double rsq = p*p + dsq, r0 = std::sqrt(p0*p0 + dsq), r1 = std::sqrt(p1*p1 + dsq);
        double f2 = std::log((r1+l1)/(r0+l0)); // Consider using atanh
        assert(!approxZero(r0+l0) && !approxZero(r1+l1));

        scaRad += phat.dot(uhat) * (
            p*f2 - d * (std::atan2(p*l1, rsq+d*r1) - std::atan2(p*l0, rsq+d*r0)));

        vecRad += uhat * (rsq*f2 + l1*r1 - l0*r0);
    }

    scaRad /= 2.0*area;
    vecRad /= 4.0*area;

    return std::make_tuple(scaRad, vecRad, R);
}

void Mesh::Triangle::buildSelfIntegrated() {
    const auto& [V0, V1, V2] = getVerts();

    double a00 = V0.dot(V0), a01 = V0.dot(V1), a02 = V0.dot(V2);
    double a11 = V1.dot(V1), a12 = V1.dot(V2), a22 = V2.dot(V2);

    double l0 = (V1-V2).norm(), l1 = (V2-V0).norm(), l2 = (V0-V1).norm();
    
    auto logN = [](double l0, double l1, double l2) {
        double lsum = l0+l1, ldiff = l2-l0;
        return std::log((lsum*lsum - l2*l2) / (l1*l1 - ldiff*ldiff));
    };

    double log0 = logN(l0, l1, l2), log1 = logN(l1, l2, l0), log2 = logN(l2, l0, l1);
    double log3 = logN(l0, l2, l1), log4 = logN(l2, l1, l0), log5 = logN(l1, l0, l2);

    // Integral of \lambda_i * \lambda_i' / r
    auto f0 = [](double l0, double l1, double l2, double log0, double log1, double log2) {
        double l0sq = l0*l0, l1sq = l1*l1, l2sq = l2*l2;
        return log0 / (20.0*l0)
            + log1 * (l0sq + 5.0*l1sq - l2sq) / (120.0*l1sq*l1)
            + log2 * (l0sq - l1sq + 5.0*l2sq) / (120.0*l2sq*l2)
            + (l2-l0) / (60.0*l1sq) + (l1-l0) / (60.0*l2sq);
    };

    // Integral of \lambda_i * \lambda_j' / r
    auto f1 = [](double l0, double l1, double l2, double log0, double log1, double log2) {
        double l0sq = l0*l0, l1sq = l1*l1, l2sq = l2*l2;
        return log2 / (40.0*l2)
            + log0 * (3.0*l0sq + l1sq - l2sq) / (80.0*l0sq*l0)
            + log1 * (l0sq + 3.0*l1sq - l2sq) / (80.0*l1sq*l1)
            + (l2-l1) / (40.0*l0sq) + (l2-l0) / (40.0*l1sq);
    };

    // Integral of \lambda_i / r or \lambda_i' / r
    auto f2 = [](double l0, double l1, double l2, double log0, double log1, double log2) {
        double l0sq = l0*l0, l1sq = l1*l1, l2sq = l2*l2;
        return log0 / (8.0*l0)
            + log1 * (l0sq + 5.0*l1sq - l2sq) / (48.0*l1sq*l1)
            + log2 * (l0sq - l1sq + 5.0*l2sq) / (48.0*l2sq*l2)
            + (l2-l0) / (24.0*l1sq) + (l1-l0) / (24.0*l2sq);
    };

    selfInts[0] 
        = f0(l0, l1, l2, log0, log1, log2) * (a00 - 2.0*a01 + a11)
        + f0(l1, l2, l0, log1, log2, log0) * (a00 - 2.0*a02 + a22)
        + f1(l0, l1, l2, log0, log1, log2) * (a00 - a02 - a01 + a12)
        + f1(l0, l2, l1, log3, log4, log5) * (a00 - a01 - a02 + a12);

    selfInts[1] = f2(l0, l1, l2, log0, log1, log2); 

    selfInts[2] = f2(l1, l2, l0, log1, log2, log0);

    selfInts[3] = (log0 / l0 + log1 / l1 + log2 / l2) / 3.0;
}

/*
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
}*/

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

