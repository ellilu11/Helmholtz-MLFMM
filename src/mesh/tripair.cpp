#include "trimoments.h"

/* TriPair(iTris)
 * Construct triangle pair from global indices of triangles
 * iTris : global indices of triangles (ordered)
 */
Mesh::TriPair::TriPair(pair2i iTris)
    : iTris(iTris)
{
    auto [iTri0, iTri1] = iTris;
    assert(iTri0 <= iTri1);
    iPair = iTri0 + iTri1*(iTri1+1)/2; // index in glTriPairs for this pair

    buildNumCommon();

    if (nCommon >= nCommonThres) {
        buildIntegratedInvR();
        if (config.alpha != 1.0) buildIntegratedInvRcubed();
    }
}

/* buildNumCommon()
 * Compute number of common vertices between triangle pair
 */
void Mesh::TriPair::buildNumCommon() {
    if (iTris.first == iTris.second) {
        nCommon = 3;
        dist = 0.0;
        return;
    }

    const auto [tri0, tri1] = getTriPair();
    nCommon = 0;
    dist = (tri0.center - tri1.center).norm();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (tri0.iVerts[i] == tri1.iVerts[j])
                ++nCommon;
    assert(nCommon < 3);
}

/* buildIntegratedInvR()
 * Build integrated 1/R and its symmetric case for triangle pair
 */
void Mesh::TriPair::buildIntegratedInvR() {
    const auto [obsTri, srcTri] = getTriPair();

    intInvR.reserve(obsTri.triQuads.size());
    for (const auto& [obs, weight] : obsTri.triQuads)
        intInvR.emplace_back(srcTri.getIntegratedInvR(obs));

    intInvR2.reserve(srcTri.triQuads.size());
    for (const auto& [src, weight] : srcTri.triQuads)
        intInvR2.emplace_back(obsTri.getIntegratedInvR(src));
}

/* buildIntegratedInvRcubed()
 * Build integrated 1/R^3 and its symmetric case for triangle pair
 */
void Mesh::TriPair::buildIntegratedInvRcubed() {
    const auto [obsTri, srcTri] = getTriPair();

    intInvRcubed.reserve(obsTri.triQuads.size());
    for (const auto& [obs, weight] : obsTri.triQuads)
        intInvRcubed.emplace_back(srcTri.getIntegratedInvRcubed(obs));

    intInvRcubed2.reserve(srcTri.triQuads.size());
    for (const auto& [src, weight] : srcTri.triQuads)
        intInvRcubed2.emplace_back(obsTri.getIntegratedInvRcubed(src));
}

Mesh::TriMoments::TriMoments(size_t nPair)
{
    nPair = glTriPairs.size();
    assert(nPair > 0);

    if (config.ie == IE::EFIE || config.ie == IE::CFIE)
        buildMomentsEFIE();

    if (config.ie == IE::MFIE || config.ie == IE::CFIE)
        buildMomentsMFIE();
}

/* buildMomentsEFIE()
 * Build EFIE moments for triangle pair
 */
void Mesh::TriMoments::buildMomentsEFIE() {
    double k = config.k;

    momentsEFIE.resize(nPair);
    for (size_t iPair = 0; iPair < nPair; ++iPair) {
        auto& [m00, m10, m01, m11] = momentsEFIE[iPair];

        const auto& triPair = glTriPairs[iPair];
        const auto& [obsTri, srcTri] = triPair.getTriPair();

    for (const auto& [obs, obsWeight] : obsTri.quads) {
        for (const auto& [src, srcWeight] : srcTri.quads) {
            double r = (obs-src).norm();

                cmplx G = obsWeight*srcWeight / (4.0*PI); // apply 1/(4pi) factor
                if (triPair.nCommon >= nCommonThres)
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
}

/* buildMomentsMFIE()
 * Build N-MFIE moments for triangle pair
 */
void Mesh::TriMoments::buildMomentsMFIE() {
    double k = config.k, k2 = k*k;

    momentsMFIE.resize(nPair);
    momentsMFIE2.resize(nPair);
    for (size_t iPair = 0; iPair < nPair; ++iPair) {
        auto& [m000, m001, m10, m01, m11] = momentsMFIE[iPair];
        auto& [n000, n001, n10, n01, n11] = momentsMFIE2[iPair];

        const auto& triPair = glTriPairs[iPair];
        const auto& [obsTri, srcTri] = triPair.getTriPair();

        vec3d obsNhat = obsTri.nhat, srcNhat = srcTri.nhat;

    for (const auto& [obs, obsWeight] : obsTri.quads) {
        // For common triangles, use analytic integration of -1/2 J term
        if (nCommon == 3) continue;

        for (const auto& [src, srcWeight] : srcTri.quads) {
            const vec3d& R = obs-src;
            double r = R.norm(), r2 = r*r, r3 = r*r2;
            assert(!Math::fzero(r));

                vec3cd gradG = obsWeight*srcWeight * R / (4.0*PI*r3); // apply 1/(4pi) factor

                if (triPair.nCommon >= nCommonThres)
                    gradG = gradG * ((-1.0+iu*k*r)*exp(iu*k*r) + 0.5*k2*r2 + 1.0); // double check signs
                else
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

/* buildIntegratedInvRcubed()
 * Build integrated 1/R^3 and its symmetric case for triangle pair
 */
void Mesh::TriPair::buildIntegratedInvRcubed() {
    const auto [obsTri, srcTri] = getTriPair();

    integratedInvRcubed.reserve(obsTri.quads.size());
    for (const auto& [obs, weight] : obsTri.quads)
        integratedInvRcubed.emplace_back(srcTri.getIntegratedInvRcubed(obs));

    integratedInvRcubed2.reserve(srcTri.quads.size());
    for (const auto& [src, weight] : srcTri.quads)
        integratedInvRcubed2.emplace_back(obsTri.getIntegratedInvRcubed(src));
}

/*
void Mesh::TriPair::buildMomentsMFIE_T() {
    momentsMFIE_T = { vec3cd::Zero(), vec3cd::Zero(), vec3cd::Zero(), 0.0 };
    auto& [m00, m10, m01, m11] = momentsMFIE_T;
    const auto& [obsTri, srcTri] = getTriPair();
    double k = config.k, k2 = k*k;
    vec3d nhat = obsTri.nhat;

    for (const auto& [obs, obsWeight] : obsTri.quads) {
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

        for (const auto& [src, srcWeight] : srcTri.quads) {
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