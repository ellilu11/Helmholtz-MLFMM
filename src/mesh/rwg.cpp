#include "rwg.h"

Mesh::RWG::RWG(
    std::shared_ptr<Exc::PlaneWave> Einc, size_t iSrc, const vec4i& idx4)
    : Source(std::move(Einc), iSrc),
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

    buildVoltage();
}

void Mesh::RWG::buildVoltage() {
    auto [inc, incNormal] = getIntegratedPlaneWave(Einc->wavevec);

    // inc . (nhat x pol) = pol . (inc x nhat) = -pol . (nhat x inc) = -pol . incNormal
    voltage = -Einc->amplitude * Einc->pol.dot(
        config.alpha * inc - config.beta * incNormal);
}

/* getIntegratedPlaneWave(kvec, doNumeric)
 * Return integral of exp(ik dot r'} * f(r') dr' at this RWG
 */
std::pair<vec3cd,vec3cd> 
Mesh::RWG::getIntegratedPlaneWave(const vec3d& kvec, bool doNumeric) const {
    vec3cd rad = vec3cd::Zero();
    vec3cd radNormal = vec3cd::Zero();
    int iTri = 0;

    if (doNumeric) {
        for (const auto& [tri, vert] : getTrisAndVerts()) {
            for (const auto& [node, weight] : tri.triQuads) {
                vec3cd nodeRad = 
                    weight * exp(iu*kvec.dot(node)) * (node - vert) * Math::sign(iTri);
                rad += nodeRad;
                radNormal += tri.nhat.cross(nodeRad);
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
        radNormal += tri.nhat.cross(triRad).conjugate(); // double check sign and conjugation!
    }

    assert(!std::isnan(rad.norm()));
    return { leng * rad, leng * radNormal };
}

/* getIntegratedEFIE(src)
 * Return the radiated field due to src tested with this RWG
 */
cmplx Mesh::RWG::getIntegratedEFIE(const std::shared_ptr<Source> src) const {
    if (config.alpha == 0) return 0.0;
    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    double k2 = config.k * config.k;
    
    cmplx intRad = 0.0;

    int iObs = 0;
    for (const auto& [obsTri, vobs] : getTrisAndVerts() ) {

        int iSrc = 0;
        for (const auto& [srcTri, vsrc] : srcRWG->getTrisAndVerts() ) {
            const TriPair& triPair = glTriPairs.at(std::minmax(obsTri.iTri, srcTri.iTri));

            const auto& [m00, m10, m01, m11] = triPair.momentsEFIE;
            vec3d v0 = (obsTri.iTri <= srcTri.iTri) ? vobs : vsrc;
            vec3d v1 = (obsTri.iTri <= srcTri.iTri) ? vsrc : vobs;
            cmplx pairRad = m11 - v1.dot(m10) - v0.dot(m01) + (v0.dot(v1) - 4.0/k2)*m00;

            // Integrate 1/R term (analytically)
            // Average obs-src and src-obs to preserve symmetry
            if (triPair.nCommon >= nCommonThres)
                pairRad += (
                    obsTri.getSingularEFIE(srcTri, triPair, vobs, vsrc) +
                    srcTri.getSingularEFIE(obsTri, triPair, vsrc, vobs)) / 2.0;

            /* For common triangles, integrate 1/R term (analytically)
            if (nCommon == 3)
                pairRad += obsTri.getDoubleSelfIntegratedInvR(vobs, vsrc);
            */

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
    if (config.alpha == 1) return 0.0;
    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    double k = config.k, k2 = k*k;

    cmplx intRad = 0.0;
    int iObs = 0;
    for (const auto& [obsTri, vobs] : getTrisAndVerts()) {

        int iSrc = 0;
        for (const auto& [srcTri, vsrc] : srcRWG->getTrisAndVerts()) {
            if (obsTri.iTri == srcTri.iTri) {
                ++iSrc;
                continue; // defer to getIntegratedMass()
            }

            const TriPair& triPair = glTriPairs.at(std::minmax(obsTri.iTri, srcTri.iTri));
            
            const auto& [m000, m001, m10, m01, m11] =
                (obsTri.iTri <= srcTri.iTri) ? triPair.momentsMFIE : triPair.momentsMFIE2;
            cmplx pairRad = m11 - vsrc.dot(m10) - vobs.dot(m01) + 
                (vobs.dot(vsrc))*m000 + obsTri.nhat.dot(vsrc)*vobs.dot(m001);

            // Integrate 1/R term (analytically)
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
}

double Mesh::RWG::getIntegratedMass(const std::shared_ptr<Source> src) const {
    if (config.alpha == 1) return 0.0;
    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    double k = config.k, k2 = k*k;

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
            for (const auto& [node, weight] : obsTri.triQuads)
                pairMass += (node-vobs).dot(node-vsrc) * weight;
            */

            // Analytic integration
            auto [X0, X1, X2] = obsTri.getVerts();
            double pairMass =
                (X0.squaredNorm() + X1.squaredNorm() + X2.squaredNorm() +
                    X0.dot(X1) + X0.dot(X2) + X1.dot(X2)) / 12.0
                - (vobs + vsrc).dot(obsTri.center) / 2.0
                + vobs.dot(vsrc) / 2.0;

            mass += pairMass / (2.0 * obsTri.area) * Math::sign(iObs) * Math::sign(iSrc);
            ++iSrc;
        }

        ++iObs;
    }

    assert(!std::isnan(mass));
    return -0.5 * leng * srcRWG->leng * mass; // factor of -1/2 = -Omega/(4pi)
}

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
            for (const auto& [obs, obsWeight] : obsTri.triQuads) {
                // For common triangles, use -1/2 J term (numerically)
                if (triPair.nCommon == 3) continue;

                for (const auto& [src, srcWeight] : srcTri.triQuads) {
                    vec3d R = obs-src;
                    double r = R.norm(), r2 = r*r, r3 = r*r2;
                    assert(!Math::fzero(r));

                    vec3cd gradG = R / r3;
                    if (triPair.nCommon == 2)
                        gradG = gradG * ((-1.0+iu*k*r)*exp(iu*k*r)+1.0+0.5*k*k*r2);
                    else if (triPair.nCommon < 2)
                        gradG = gradG * (-1.0+iu*k*r)*exp(iu*k*r);

                    // minus sign from flipping J x gradG to gradG x J
                    pairRad -= (obs-vobs).dot(nhat.cross(gradG.cross(src-vsrc)))
                        * obsWeight * srcWeight;
                }
            }

            // For edge adjacent triangles, integrate 1/R term (analytically)
            if (triPair.nCommon == 2)
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