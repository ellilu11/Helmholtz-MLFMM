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

    /* Update triToRWG
    for (int i = 0; i < 2; ++i) {
        int iTri = iTris[i];
        auto& triToRWG = triToRWGs[iTri];
        triToRWG.iRWGs.push_back(iSrc);

        assert(iTri == idx4[2] || iTri == idx4[3]);
        triToRWG.isMinus.push_back(iTri == idx4[3]);

        assert(triToRWG.iRWGs.size() <= 3 && triToRWG.isMinus.size() <= 3);
    }*/
}

/* getIntegratedPlaneWave(kvec, doNumeric)
 * Return integral of exp(ik dot r'} * f(r') dr' at this RWG
 */
vec3cd Mesh::RWG::getIntegratedPlaneWave(const vec3d& kvec, bool doNumeric) const {
    using namespace Math;

    auto Xnc = getVertsNC();
    vec3cd rad = vec3cd::Zero();
    int iTri = 0;

    if (doNumeric) {
        for (const auto& tri : getTris()) {
            for (const auto& [node, weight] : tri.triQuads)
                rad += weight * exp(iu*kvec.dot(node))
                        * (node - Xnc[iTri])
                        * sign(iTri);
            ++iTri;
        }

        return leng * rad;
    }

    for (const auto& tri : getTris()) {
        const vec3d& X0 = tri.getVerts()[0];
        const auto& [scaRad, vecRad] = tri.getIntegratedPlaneWave(kvec);
        rad += exp(iu*kvec.dot(X0)) // TODO: Address warning
                * (scaRad * (X0 - Xnc[iTri]) + vecRad) 
                * sign(iTri++);
    }

    assert(!std::isnan(rad.norm()));
    return leng * rad;
}

/* getIntegratedEFIE(src)
 * Return the electric field due to src tested with this RWG
 */
cmplx Mesh::RWG::getIntegratedEFIE(const std::shared_ptr<Source> src) const {
    if (config.alpha == 0) return 0.0;
    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    cmplx intRad = 0.0;

    int iObs = 0;
    for (const auto& [obsTri, vobs] : getTrisAndVerts() ) {

        int iSrc = 0;
        for (const auto& [srcTri, vsrc] : srcRWG->getTrisAndVerts() ) {
            vec3d v0 = (obsTri.iTri <= srcTri.iTri) ? vobs : vsrc;
            vec3d v1 = (obsTri.iTri <= srcTri.iTri) ? vsrc : vobs;

            const TriPair& triPair = glTriPairs.at(std::minmax(obsTri.iTri, srcTri.iTri));
            const auto& [m00, m10, m01, m11] = triPair.momentsEFIE;

            cmplx pairRad =
                m11 - v1.dot(m10) - v0.dot(m01) + (v0.dot(v1) - 4.0/(k*k))*m00;

            // For edge adjacent triangles, integrate 1/R term (analytically)
            // Average obs-src and src-obs to preserve symmetry
            if (triPair.nCommon >= 2)
                pairRad += (
                    obsTri.getDoubleIntegratedSingularEFIE(srcTri, triPair, vobs, vsrc) +
                    srcTri.getDoubleIntegratedSingularEFIE(obsTri, triPair, vsrc, vobs)) / 2.0;

            /* For common triangles, integrate 1/R term (analytically)
            if (nCommon == 3)
                pairRad += obsTri.getDoubleSelfIntegratedInvR(vobs, vsrc);
            */
           
            intRad += pairRad * Math::sign(iObs) * Math::sign(iSrc);
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
    cmplx intRad = 0.0;

    /*
    int iObsTri = 0;
    for (const auto& obsTri : getTris()) {
        const auto& obsNC = getVertsNC()[iObsTri];

        int iSrcTri = 0;
        for (const auto& srcTri : srcRWG->getTris()) {
            int nCommon = obsTri.getNumCommonVerts(srcTri);
            const vec3d& srcNC = srcRWG->getVertsNC()[iSrcTri];

            cmplx pairRad = 0.0;
            for (const auto& [obs, obsWeight] : obsTri.triQuads) {
                for (const auto& [src, srcWeight] : srcTri.triQuads) {
                    const vec3d& rvec = obs-src;
                    double r = rvec.norm(), r2 = r*r, r3 = r*r2;

                    vec3d gradG = rvec / r3;
                    if (nCommon == 2) 
                        gradG = gradG * ((-1.0+iu*k*r)*exp(iu*k*r)+1.0+0.5*k*k*r2); // double check signs
                    else if (nCommon < 2)
                        gradG = gradG * (-1.0+iu*k*r)*exp(iu*k*r);

                    pairRad -= (obs-obsNC).dot(gradG.cross(src-srcNC))
                                * obsWeight * srcWeight;
                }

                if (nCommon == 3) {
                    pairRad -= 0.5 * (obs-obsNC).dot(obsTri.nhat.cross(obs-obsNC))
                        * obsWeight;
                }
            }

            if (nCommon == 2)
                pairRad += obsTri.getDoubleIntegratedSingularMFIE(srcTri, obsNC, srcNC);

            intRad += pairRad * Math::sign(iSrcTri) * Math::sign(iObsTri);
            ++iSrcTri;
        }

        ++iObsTri;
    }*/

    assert(!std::isnan(intRad.real()) && !std::isnan(intRad.imag()));
    return leng * srcRWG->leng * intRad;
}
