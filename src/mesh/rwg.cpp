#include "rwg.h"

Mesh::RWG::RWG(const vec4i& idx4, size_t iSrc)
    : Source(iSrc),
    iTris({ idx4[2], idx4[3] }),
    iVertsC({ idx4[0], idx4[1] })
{
    auto tris = getTris();
    auto verts = getVertsC();

    leng = (verts[0]-verts[1]).norm();

    // Find indices of non-common vertices
    for (int i = 0; i < 2; ++i)
        for (const auto& iVert : tris[i].iVerts)
            if (iVert != iVertsC[0] && iVert != iVertsC[1])
                iVertsNC[i] = iVert;
}

/* getVoltage()
 * Return voltage of this RWG due to superposition of incident plane waves
 */
cmplx Mesh::RWG::getVoltage() {
    using namespace Exct;

    cmplx voltage = 0.0;
    for (const auto& Einc : Eincs) {
        auto [inc, incNormal] = getIntegratedPlaneWave(Einc->wavevec);

        // (1-alpha) * eta * H_inc = (1-alpha) * khat x E_inc
        vec3d polE = config.alpha * Einc->pol;
        vec3d polH = (1.0 - config.alpha) * Einc->wavehat.cross(Einc->pol);

        // inc . (nhat x polH) = polH . (inc x nhat) = -polH . (nhat x inc) = -polH . incNormal
        voltage -= Einc->amplitude * (polE.dot(inc) - polH.dot(incNormal));
    }

    return voltage;
}

/* getIntegratedPlaneWave(kvec, doNumeric)
 * Return integral of exp(i kvec . r'} * f(r') dr' at this RWG
 */
std::pair<vec3cd,vec3cd> 
Mesh::RWG::getIntegratedPlaneWave(const vec3d& kvec, bool doNumeric) const {
    vec3cd rad = vec3cd::Zero();
    vec3cd radNormal = vec3cd::Zero();
    int iTri = 0;

    if (doNumeric) {
        for (const auto& [tri, vert] : getTrisAndVerts()) {
            for (const auto& [node, weight] : tri.quads) {
                vec3cd nodeRad = 
                    weight * exp(iu*kvec.dot(node)) * (node - vert) * Math::sign(iTri);
                rad += nodeRad;
                if (config.ie != IE::EFIE) 
                    radNormal += tri.nhat.cross(nodeRad).conjugate();
            }
            ++iTri;
        }
        return { leng * rad, leng * radNormal };
    }

    for (const auto& [tri, vert] : getTrisAndVerts()) {
        vec3d X0 = tri.getVerts()[0];
        const auto& [scaRad, vecRad] = tri.getIntegratedPlaneWave(kvec);

        vec3cd triRad = 
            exp(iu*kvec.dot(X0)) * (scaRad * (X0 - vert) + vecRad) * Math::sign(iTri++);

        rad += triRad;
        if (config.ie != IE::EFIE)
            radNormal += tri.nhat.cross(triRad).conjugate(); // double check sign and conjugation!
    }

    assert(!std::isnan(rad.norm()));
    return { leng * rad, leng * radNormal };
}

/* getIntegratedEFIE(src)
 * Return the electric field due to src tested with this RWG
 */
cmplx Mesh::RWG::getIntegratedEFIE(const std::shared_ptr<Source> src) const {
    if (config.ie == IE::MFIE) return 0.0;
    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    double k2 = config.k * config.k;
    
    cmplx intRad = 0.0;

    int iObs = 0;
    for (const auto& [obsTri, vobs] : getTrisAndVerts() ) {

        int iSrc = 0;
        for (const auto& [srcTri, vsrc] : srcRWG->getTrisAndVerts() ) {
            size_t iPair = glPairsToIdx.at(std::minmax(obsTri.iTri, srcTri.iTri));

            const auto& [m00, m10, m01, m11] = glTriPairs.momentsEFIE[iPair];
            vec3d v0 = (obsTri.iTri <= srcTri.iTri) ? vobs : vsrc;
            vec3d v1 = (obsTri.iTri <= srcTri.iTri) ? vsrc : vobs;
            cmplx pairRad = m11 - v1.dot(m10) - v0.dot(m01) + (v0.dot(v1) - 4.0/k2)*m00;

            // Integrate singular term (numeric-analytic)
            // Average obs-src and src-obs to preserve symmetry
            if (glTriPairs.nCommons[iPair] >= nCommonThres)
                pairRad += (
                    obsTri.getSingularEFIE(srcTri, vobs, vsrc, iPair) +
                    srcTri.getSingularEFIE(obsTri, vsrc, vobs, iPair)) / 2.0;

            // For common triangles, integrate 1/R term (full-analytic)
            //if (triPair.nCommon == 3)
            //    pairRad += obsTri.getDoubleSelfIntegratedInvR(vobs, vsrc);

            intRad += pairRad * Math::sign(iSrc) * Math::sign(iObs);
            ++iSrc;
        }

        ++iObs;
    }

    assert(!std::isnan(intRad.real()) && !std::isnan(intRad.imag()));
    return leng * srcRWG->leng * intRad;
}

/* getIntegratedMFIE(src)
 * Return the magnetic field due to src tested with this RWG
 */
cmplx Mesh::RWG::getIntegratedMFIE(const std::shared_ptr<Source> src) const {
    if (config.ie == IE::EFIE) return 0.0;
    const auto srcRWG = dynamic_pointer_cast<RWG>(src);

    cmplx intRad = 0.0;
    int iObs = 0;
    for (const auto& [obsTri, vobs] : getTrisAndVerts()) {

        int iSrc = 0;
        for (const auto& [srcTri, vsrc] : srcRWG->getTrisAndVerts()) {
            if (obsTri.iTri == srcTri.iTri) {
                ++iSrc;
                continue; // defer to getIntegratedMass()
            }

            size_t iPair = glPairsToIdx.at(std::minmax(obsTri.iTri, srcTri.iTri));
            
            const auto& [m000, m001, m10, m01, m11] =
                (obsTri.iTri <= srcTri.iTri) ? 
                    glTriPairs.momentsMFIE[iPair] : glTriPairs.momentsMFIE2[iPair];
            cmplx pairRad = m11 - vsrc.dot(m10) - vobs.dot(m01) + 
                (vobs.dot(vsrc))*m000 + obsTri.nhat.dot(vsrc)*vobs.dot(m001);

            // Integrate singular terms (numeric-analytic)
            if (glTriPairs.nCommons[iPair] >= nCommonThres)
                pairRad += // double check sign!
                    obsTri.getSingularMFIE(srcTri, vobs, vsrc, iPair);

            intRad += pairRad * Math::sign(iObs) * Math::sign(iSrc);
            ++iSrc;
        }

        ++iObs;
    }

    assert(!std::isnan(intRad.real()) && !std::isnan(intRad.imag()));
    return leng * srcRWG->leng * intRad;
}

/* getIntegratedMass(src)
 * Return the mass matrix entry for this RWG and src
 */
double Mesh::RWG::getIntegratedMass(const std::shared_ptr<Source> src) const {
    if (config.ie == IE::EFIE) return 0.0;
    const auto srcRWG = dynamic_pointer_cast<RWG>(src);

    double mass = 0.0;
    int iObs = 0;
    for (const auto& [obsTri, vobs] : getTrisAndVerts()) {

        int iSrc = 0;
        for (const auto& [srcTri, vsrc] : srcRWG->getTrisAndVerts()) {
            if (obsTri.iTri != srcTri.iTri) {
                ++iSrc;
                continue; // defer to getIntegratedMFIE()
            }

            /* Numeric integration
            double pairMass = 0.0;
            for (const auto& [node, weight] : obsTri.quads)
                pairMass += (node-vobs).dot(node-vsrc) * weight;
            */

            // Analytic integration
            auto [X0, X1, X2] = obsTri.getVerts();
            double pairMass =
                (X0.squaredNorm() + X1.squaredNorm() + X2.squaredNorm() +
                    X0.dot(X1) + X0.dot(X2) + X1.dot(X2)) / 12.0
                - (vobs + vsrc).dot(obsTri.center) / 2.0
                + vobs.dot(vsrc) / 2.0;

            mass += pairMass / (2.0 * obsTri.area) 
                * Math::sign(iObs) * Math::sign(iSrc);
            ++iSrc;
        }

        ++iObs;
    }

    assert(!std::isnan(mass));
    // Puzzle: Why plus sign here?
    return 0.5 * leng * srcRWG->leng * mass; // factor of 1/2 = Omega/(4pi)
}

/* Using inc . (nhat x polH) directly
void Mesh::RWG::buildVoltage() {
    vec3d kvec = Einc->wavevec;
    vec3d polE = config.alpha * Einc->pol;
    vec3d polH = (1.0 - config.alpha) * Einc->pol; // Einc->wavehat.cross(Einc->pol);

    vec3cd inc = vec3cd::Zero();
    cmplx voltH = 0.0;
    int iTri = 0;
    for (const auto& [tri, vert] : getTrisAndVerts()) {
        vec3d X0 = tri.getVerts()[0];
        const auto& [scaRad, vecRad] = tri.getIntegratedPlaneWave(kvec);

        vec3cd triRad =
            exp(iu*kvec.dot(X0)) * (scaRad * (X0 - vert) + vecRad) * Math::sign(iTri++);

        inc += triRad;
        voltH += tri.nhat.cross(polH).dot(triRad); // Hermitian dot!
    }

    inc *= leng;
    voltH *= leng;

    voltage = -Einc->amplitude * (polE.dot(inc) + voltH);
}
*/

/* Debug version, no precomputed moments
cmplx Mesh::RWG::getIntegratedMFIE(const std::shared_ptr<Source> src) const {
    if (config.alpha == 1.0) return 0.0;
    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    double k = config.k, k2 = k*k;

    cmplx intRad = 0.0;
    int iObs = 0;
    for (const auto& [obsTri, vobs] : getTrisAndVerts()) {
        vec3d nhat = obsTri.nhat;

        int iSrc = 0;
        for (const auto& [srcTri, vsrc] : srcRWG->getTrisAndVerts()) {
            const TriPair& triPair = glTriPairs.at(std::minmax(obsTri.iTri, srcTri.iTri));

            cmplx pairRad = 0.0;
            for (const auto& [obs, obsWeight] : obsTri.quads) {
                // For common triangles, use -1/2 J term (numerically)
                if (triPair.nCommon == 3) continue;

                for (const auto& [src, srcWeight] : srcTri.quads) {
                    vec3d R = obs-src;
                    double r = R.norm(), r2 = r*r, r3 = r*r2;
                    assert(!lzero(r));

                    vec3cd gradG = obsWeight*srcWeight * R / (4.0*PI*r3);
                    if (triPair.nCommon >= nCommonThres)
                        gradG = gradG * ((-1.0+iu*k*r)*exp(iu*k*r)+1.0+0.5*k*k*r2);
                    else
                        gradG = gradG * (-1.0+iu*k*r)*exp(iu*k*r);

                    // minus sign from flipping J x gradG to gradG x J
                    pairRad -= (obs-vobs).dot(nhat.cross(gradG.cross(src-vsrc)));
                }
            }

            // For edge adjacent triangles, integrate 1/R term (analytically)
            if (triPair.nCommon >= nCommonThres)
                pairRad -= // minus sign since 1.0+0.5*k*k*r2 was added to gradG in MFIE moments
                    obsTri.getSingularMFIE(srcTri, triPair, vobs, vsrc);

            intRad += pairRad * Math::sign(iObs) * Math::sign(iSrc);
            ++iSrc;
        }

        ++iObs;
    }

    assert(!std::isnan(intRad.real()) && !std::isnan(intRad.imag()));
    return leng * srcRWG->leng * intRad;
}*/