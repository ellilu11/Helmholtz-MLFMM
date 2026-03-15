#include "tripair.h"

/* TriPair(iTris)
 * Construct triangle pair from global indices of triangles
 * iTris : global indices of triangles (ordered)
 */
Mesh::TriPair::TriPair(pair2i iTris)
    : iTris(iTris)
{
    assert(iTris.first <= iTris.second);
    buildNumCommon();
    buildMomentsEFIE();
    if (nCommon >= 2) buildIntegratedInvR();
}

/* buildNumCommon()
 * Compute number of common vertices between triangle pair
 */
void Mesh::TriPair::buildNumCommon() {
    if (iTris.first == iTris.second) {
        nCommon = 3;
        return;
    }

    const auto [tri0, tri1] = getTriPair();
    nCommon = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (tri0.iVerts[i] == tri1.iVerts[j])
                ++nCommon;
    assert(nCommon < 3);
}

/* buildMomentsEFIE()
 * Build EFIE moments for triangle pair
 */
void Mesh::TriPair::buildMomentsEFIE() {
    double k = config.k;
    const auto& [obsTri, srcTri] = getTriPair();

    momentsEFIE = { 0.0, vec3cd::Zero(), vec3cd::Zero(), 0.0 };
    for (const auto& [obs, obsWeight] : obsTri.quads) {
        for (const auto& [src, srcWeight] : srcTri.quads) {
            double r = (obs-src).norm();

            cmplx G = obsWeight * srcWeight;
            if (nCommon >= 2) 
                G *= (Math::fzero(r) ? iu*k : (exp(iu*k*r)-1.0) / r);
            else {
                assert(!Math::fzero(r));
                G *= exp(iu*k*r) / r;
            }

            auto& [m00, m10, m01, m11] = momentsEFIE;
            m00 += G;
            m10 += obs * G;
            m01 += src * G;
            m11 += obs.dot(src) * G;
        }
    }
}

/* buildIntegratedInvR()
 * Build integrated 1/R and its symmetric case for triangle pair
 */
void Mesh::TriPair::buildIntegratedInvR() {
    const auto [obsTri, srcTri] = getTriPair();

    integratedInvR.reserve(obsTri.quads.size());
    for (const auto& [obs, weight] : obsTri.quads)
        integratedInvR.emplace_back(srcTri.getIntegratedInvR(obs));

    integratedInvR2.reserve(srcTri.quads.size());
    for (const auto& [src, weight] : srcTri.quads)
        integratedInvR2.emplace_back(obsTri.getIntegratedInvR(src));
}

/*
void Mesh::TriPair::buildMomentsInvR() {
    const auto& [obsTri, srcTri] = getTriPair();

    momentsInvR = { 0.0, vec3cd::Zero(), 0.0 };
    for (const auto& [obs, obsWeight] : obsTri.quads) {
        const auto& [scaRad, vecRad] = srcTri.getIntegratedInvR(obs);

        auto& [m0, m1, m11] = momentsInvR;
        m0 += scaRad * obsWeight;
        m1 += vecRad * obsWeight;
        m11 += obs.dot(vecRad) * obsWeight;
    }

    // TODO: Symmetric case
}
*/