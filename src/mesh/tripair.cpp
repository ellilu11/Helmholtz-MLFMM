#include "tripair.h"

Mesh::TriPair::TriPair(pair2i iTris)
    : iTris(iTris)
{
    assert(iTris.first <= iTris.second);
    buildNumCommon();

    if (config.alpha != 0.0) buildMomentsEFIE();
    if (config.alpha != 1.0) buildMomentsMFIE();
    if (nCommon >= 2) {
        buildIntegratedInvR();
        if (config.alpha != 1.0) buildIntegratedInvRcubed();
    }
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
    momentsEFIE = { 0.0, vec3cd::Zero(), vec3cd::Zero(), 0.0 };
    const auto& [obsTri, srcTri] = getTriPair();

    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        for (const auto& [src, srcWeight] : srcTri.triQuads) {
            double r = (obs-src).norm();

            cmplx G = obsWeight * srcWeight; // absorb quad weights into G
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

void Mesh::TriPair::buildMomentsMFIE() {
    momentsMFIE = { vec3cd::Zero(), vec3cd::Zero(), vec3cd::Zero(), 0.0 };
    if (nCommon == 3) return; // principal value of self-interaction is zero
    const auto& [obsTri, srcTri] = getTriPair();

    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        for (const auto& [src, srcWeight] : srcTri.triQuads) {
            const vec3d& R = obs-src;
            double r = R.norm(), r2 = r*r, r3 = r*r2;
            assert(!Math::fzero(r));

            vec3cd gradG = R / r3 * (-1.0+iu*k*r)*exp(iu*k*r)
                * obsWeight * srcWeight; // absorb quad weights into gradG
            /*
            if (nCommon == 2)
                gradG = gradG * ((-1.0+iu*k*r)*exp(iu*k*r)+1.0+0.5*k*k*r2); // double check signs
            else if (nCommon < 2)
                gradG = gradG * (-1.0+iu*k*r)*exp(iu*k*r);
            */

            auto& [m00, m10, m01, m11] = momentsMFIE;
            // minus signs from flipping J x gradG to gradG x J
            m00 -= gradG; 
            m10 -= obs.cross(gradG); // double check
            m01 -= (gradG.cross(src)).conjugate();
            m11 -= conj(obs.dot(gradG.cross(src)));

            auto& [m00_2, m10_2, m01_2, m11_2] = momentsMFIE2;
            m00_2 += gradG;
            m10_2 += src.cross(gradG);
            m01_2 += (gradG.cross(obs)).conjugate();
            m11_2 += conj(src.dot(gradG.cross(obs)));
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

void Mesh::TriPair::buildIntegratedInvRcubed() {
    const auto [obsTri, srcTri] = getTriPair();

    for (const auto& [obs, weight] : obsTri.triQuads)
        integratedInvRcubed.push_back(srcTri.getIntegratedInvRcubed(obs));

    for (const auto& [src, weight] : srcTri.triQuads)
        integratedInvRcubed2.push_back(obsTri.getIntegratedInvRcubed(src));
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