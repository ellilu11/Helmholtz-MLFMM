#include "tripair.h"

Mesh::TriPair::TriPair(pair2i iTris)
    : iTris(iTris)
{
    assert(iTris.first <= iTris.second);
    buildNumCommon();

    if (config.alpha != 0.0) buildMomentsEFIE();
    if (config.alpha != 1.0) buildMomentsMFIE();
    if (nCommon >= nCommonThres) {
        buildIntegratedInvR();
        if (config.alpha != 1.0) buildIntegratedInvRcubed();
    }

     //std::cout << "Built TriPair (" << iTris.first << "," << iTris.second
     //    << ") with nCommon = " << nCommon << '\n';
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
    double k = config.k;

    auto& [m00, m10, m01, m11] = momentsEFIE;
    const auto& [obsTri, srcTri] = getTriPair();

    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        for (const auto& [src, srcWeight] : srcTri.triQuads) {
            double r = (obs-src).norm();

            cmplx G = obsWeight*srcWeight / (4.0*PI); // apply 1/(4pi) factor
            if (nCommon >= nCommonThres)
                G *= (Math::fzero(r) ? iu*k : (exp(iu*k*r)-1.0) / r);
            else {
                assert(!Math::fzero(r));
                G *= exp(iu*k*r) / r;
            }

            m00 += G;
            m10 += obs * G;
            m01 += src * G;
            m11 += obs.dot(src) * G;
        }
    }
}

// Build N-MFIE moments
void Mesh::TriPair::buildMomentsMFIE() {
    double k = config.k, k2 = k*k;

    auto& [m000, m001, m10, m01, m11] = momentsMFIE;
    auto& [n000, n001, n10, n01, n11] = momentsMFIE2;

    const auto& [obsTri, srcTri] = getTriPair();
    vec3d obsNhat = obsTri.nhat, srcNhat = srcTri.nhat;

    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        // For common triangles, use analytic integration of -1/2 J term
        if (nCommon == 3) continue;

        for (const auto& [src, srcWeight] : srcTri.triQuads) {
            const vec3d& R = obs-src;
            double r = R.norm(), r2 = r*r, r3 = r*r2;
            assert(!Math::fzero(r));

            vec3cd gradG = obsWeight*srcWeight * R / (4.0*PI*r3); // apply 1/(4pi) factor

            if (nCommon >= nCommonThres)
                gradG = gradG * ((-1.0+iu*k*r)*exp(iu*k*r) + 0.5*k2*r2 + 1.0); // double check signs
            else if (nCommon < nCommonThres)
                gradG = gradG * (-1.0+iu*k*r)*exp(iu*k*r);

            // Overall minus sign from flipping J x gradG to gradG x J
            m000 -= -obsNhat.dot(gradG);
            m001 -= gradG;
            m10 -= (obs.dot(gradG) * obsNhat - obsNhat.dot(gradG) * obs);
            m01 -= obsNhat.cross(gradG.cross(src));
            m11 -= obs.dot(obsNhat.cross(gradG.cross(src)));

            n000 += -srcNhat.dot(gradG);
            n001 += gradG;
            n10 += (src.dot(gradG) * srcNhat - srcNhat.dot(gradG) * src);
            n01 += srcNhat.cross(gradG.cross(obs));
            n11 += src.dot(srcNhat.cross(gradG.cross(obs)));
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
void Mesh::TriPair::buildMomentsMFIE_T() {
    momentsMFIE_T = { vec3cd::Zero(), vec3cd::Zero(), vec3cd::Zero(), 0.0 };
    auto& [m00, m10, m01, m11] = momentsMFIE_T;
    const auto& [obsTri, srcTri] = getTriPair();
    double k = config.k, k2 = k*k;
    vec3d nhat = obsTri.nhat;

    for (const auto& [obs, obsWeight] : obsTri.triQuads) {
        // For common triangles, use -1/2 J term
        if (nCommon == 3) {
            // since numerical integration only cancels one RWG's 1/(2A) factor
            double weight = obsWeight / (4.0*obsTri.area);
            vec3d nhat_X_obs = obsTri.nhat.cross(obs);
            m00 -= obsTri.nhat * weight;
            m10 -= -nhat_X_obs * weight;
            m01 -= nhat_X_obs * weight;
            continue;
        }

        for (const auto& [src, srcWeight] : srcTri.triQuads) {
            const vec3d& R = obs-src;
            double r = R.norm(), r2 = r*r, r3 = r*r2;
            assert(!Math::fzero(r));

            vec3cd gradG = R / r3 * obsWeight*srcWeight; // T-MFIE

            if (nCommon == 2)
                gradG = gradG * ((-1.0+iu*k*r)*exp(iu*k*r) + 0.5*k2*r2 + 1.0); // double check signs
            else if (nCommon < 2)
                gradG = gradG * (-1.0+iu*k*r)*exp(iu*k*r);

            vec3cd obs_X_gradG = obs.cross(gradG).conjugate();
            vec3cd gradG_X_src = gradG.cross(src).conjugate();
            // minus signs from flipping J x gradG to gradG x J
            m00 -= gradG;
            m10 -= obs_X_gradG;
            m01 -= gradG_X_src;
            m11 -= obs.dot(gradG_X_src);
        }
    }
}*/

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