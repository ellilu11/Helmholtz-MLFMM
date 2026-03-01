#include "tripair.h"

Mesh::TriPair::TriPair(pair2i iTris)
    : iTris(iTris)
{
    assert(iTris.first <= iTris.second);
    buildNumCommon();
    buildMomentsEFIE();
    if (nCommon >= 2) buildIntegratedInvR();
}

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

void Mesh::TriPair::buildMomentsEFIE() {
    const auto& [obsTri, srcTri] = getTriPair();

    momentsEFIE = { 0.0, vec3cd::Zero(), vec3cd::Zero(), 0.0 };
    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        for (const auto& [src, srcWeight] : srcTri.triQuads) {
            double r = (obs-src).norm();

            cmplx G = obsWeight * srcWeight;
            if (nCommon >= 2) G *= (Math::fzero(r) ? iu*k : (exp(iu*k*r)-1.0) / r);
            else G *= exp(iu*k*r) / r;

            auto& [m00, m10, m01, m11] = momentsEFIE;
            m00 += G;
            m10 += obs * G;
            m01 += src * G;
            m11 += obs.dot(src) * G;
        }
    }
}

void Mesh::TriPair::buildIntegratedInvR() {
    const auto [obsTri, srcTri] = getTriPair();

    for (const auto& [obs, weight] : obsTri.triQuads)
        integratedInvR.push_back(srcTri.getIntegratedInvR(obs));

    for (const auto& [src, weight] : srcTri.triQuads)
        integratedInvR2.push_back(obsTri.getIntegratedInvR(src));
}

/*
void Mesh::TriPair::buildMomentsInvR() {
    const auto& [obsTri, srcTri] = getTriPair();

    momentsInvR = { 0.0, vec3cd::Zero(), 0.0 };
    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        const auto& [scaRad, vecRad] = srcTri.getIntegratedInvR(obs);

        auto& [m0, m1, m11] = momentsInvR;
        m0 += scaRad * obsWeight;
        m1 += vecRad * obsWeight;
        m11 += obs.dot(vecRad) * obsWeight;
    }

    // TODO: Symmetric case
}
*/