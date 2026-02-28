#include "tripair.h"

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
    const auto& [tri0, tri1] = getTriPair();

    radMoments = { 0.0, vec3d::Zero(), vec3d::Zero(), 0.0 };
    for (const auto& [obs, obsWeight] : tri0.triQuads) {
        for (const auto& [src, srcWeight] : tri1.triQuads) {
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
    const auto& [tri0, tri1] = getTriPair();

    for (const auto& [node, weight] : tri0.triQuads)
        integratedInvR.push_back(tri1.getIntegratedInvR(node));

    // TODO: Symmetric case
}

double Mesh::TriPair::getDoubleIntegratedInvR(
    const vec3d& obsNC, const vec3d& srcNC) const
{
    const auto& [obsTri, srcTri] = getTriPair();
    const vec3d& srcNCproj = srcTri.proj(srcNC);

    double rad = 0.0;
    size_t iNode = 0;
    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        const vec3d& obsProj = srcTri.proj(obs);

        const auto& [scaRad, vecRad] = integratedInvR[iNode++];
        rad += ((obs-obsNC).dot(vecRad+(obsProj-srcNCproj)*scaRad) - 4.0/(k*k)*scaRad)
            * obsWeight;
    }

    return rad;
}