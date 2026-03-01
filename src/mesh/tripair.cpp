#include "tripair.h"

Mesh::TriPair::TriPair(pair2i iTris)
    : iTris(iTris)
{
    assert(iTris.first <= iTris.second);
    buildNumCommon();
    buildRadMoments();
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

void Mesh::TriPair::buildRadMoments() {
    const auto& [obsTri, srcTri] = getTriPair();

    radMoments = { 0.0, vec3cd::Zero(), vec3cd::Zero(), 0.0 };
    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        for (const auto& [src, srcWeight] : srcTri.triQuads) {
            double r = (obs-src).norm();

            cmplx G;
            if (nCommon >= 2) G = (Math::fzero(r) ? iu*k : (exp(iu*k*r)-1.0) / r);
            else G = exp(iu*k*r) / r;
            G *= obsWeight * srcWeight;

            auto& [m00, m10, m01, m11] = radMoments;
            m00 += G;
            m10 += obs * G;
            m01 += src * G;
            m11 += obs.dot(src) * G;
        }
    }
}

void Mesh::TriPair::buildIntegratedInvR() {
    const auto& [obsTri, srcTri] = getTriPair();

    for (const auto& [obs, weight] : obsTri.triQuads)
        integratedInvR.push_back(srcTri.getIntegratedInvR(obs));

    for (const auto& [src, weight] : srcTri.triQuads)
        integratedInvR2.push_back(obsTri.getIntegratedInvR(src));
}

double Mesh::TriPair::getDoubleIntegratedInvR(const vec3d& vobs, const vec3d& vsrc, bool isSym) const
{
    const auto& [obsTri, srcTri] = getTriPair();
    const vec3d& vsrcProj = srcTri.proj(vsrc);

    double rad = 0.0;
    size_t iNode = 0;
    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        const vec3d& obsProj = srcTri.proj(obs);

        const auto& [scaRad, vecRad] = 
            (isSym ? 
            integratedInvR2[iNode] : 
            integratedInvR[iNode]);
        rad += ((obs-vobs).dot(vecRad+(obsProj-vsrcProj)*scaRad) - 4.0/(k*k)*scaRad)
            * obsWeight;

        ++iNode;
    }

    return rad;
}