#include "rwg.h"

Mesh::RWG::RWG(
    std::shared_ptr<Exc::PlaneWave> Einc, size_t iSrc, const vec4i& idx4)
    : Source(std::move(Einc), iSrc),
    iTris({ idx4[2], idx4[3] }),
    iVertsC({ idx4[0], idx4[1] })
{
    const auto& tris = getTris();
    const auto& verts = getVertsC();

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

/* getPlaneWaveIntegrated(kvec, doNumeric)
 * Return integral of exp(ik dot r'} * f(r') dr' at this RWG
 */
vec3cd Mesh::RWG::getPlaneWaveIntegrated(const vec3d& kvec, bool doNumeric) const {
    using namespace Math;

    const auto& Xnc = getVertsNC();
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
        const auto& [scaRad, vecRad] = tri.getPlaneWaveIntegrated(kvec);
        rad += exp(iu*kvec.dot(X0)) 
                * (scaRad * (X0 - Xnc[iTri]) + vecRad) 
                * sign(iTri++);
    }

    assert(!std::isnan(rad.norm()));
    return leng * rad;
}

/* getIntegratedRad(src)
 * Return the radiated field due to src tested with this RWG
 */
cmplx Mesh::RWG::getIntegratedRad(const std::shared_ptr<Source> src) const {
    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    cmplx intRad = 0.0;

    int iObsTri = 0;
    for (const auto& obsTri : getTris() ) {
        const auto& obsNC = getVertsNC()[iObsTri];

        int iSrcTri = 0;
        for (const auto& srcTri : srcRWG->getTris() ) {
            const auto 
                &srcNC = srcRWG->getVertsNC()[iSrcTri],
                &srcNCproj = srcTri.proj(srcNC);
            int nCommon = obsTri.getNumCommonVerts(srcTri);

            cmplx pairRad = 0.0;
            for (const auto& [obs, obsWeight] : obsTri.triQuads) {
                const auto& obsProj = srcTri.proj(obs);

                // Add contribution from e^{ikR)/R or (e^(ikR)-1)/R term (numerically)
                for (const auto& [src, srcWeight] : srcTri.triQuads) {
                    double r = (obs-src).norm();

                    cmplx G;
                    if (nCommon >= 2) G = (Math::fzero(r) ? iu*k : (exp(iu*k*r)-1.0) / r); 
                    else G = exp(iu*k*r) / r;

                    pairRad += ((obs-obsNC).dot(src-srcNC) - 4.0 / (k*k)) * G
                        * obsWeight * srcWeight;
                }

                // For edge adjacent triangles, add contribution from 1/R term (analytically)
                if (nCommon >= 2) {
                    const auto& [scaRad, vecRad] = srcTri.getNearIntegrated(obs);
                    pairRad +=
                        ((obs-obsNC).dot(vecRad+(obsProj-srcNCproj)*scaRad) - 4.0/(k*k)*scaRad)
                        * obsWeight;
                }
            }

            /* For common triangles, add contribution from 1/R term (analytically)
            if (nCommon == 3) {
                const auto [V0, V1, V2] = obsTri.getVerts();
                double a00 = V0.dot(V0), a01 = V0.dot(V1), a02 = V0.dot(V2); // cache?
                const auto& sumNC = srcNC + obsNC;

                pairRad += obsTri.selfInts[0]
                    + obsTri.selfInts[1] * (-2.0*a00 + 2.0*a01 + (V0-V1).dot(sumNC))
                    + obsTri.selfInts[2] * (-2.0*a00 + 2.0*a02 + (V0-V2).dot(sumNC))
                    + obsTri.selfInts[3] * (a00 - V0.dot(sumNC) + srcNC.dot(obsNC) - 4.0/(k*k));
            }
            */
            
            // Using precomputed moments
            //const auto [mm0, mm1, mm2, mm3] =
            //    Mesh::glRadMoments.at(makeUnordered(obsTri.iTri, srcTri.iTri));
            //intRad +=
            //    (mm0 - obsNC.dot(mm1) - mm2.dot(srcNC) + mm3 * (obsNC.dot(srcNC) - 4.0/(k*k)));

            intRad += pairRad * Math::sign(iSrcTri) * Math::sign(iObsTri);
            ++iSrcTri;
        }

        ++iObsTri;
    }

    assert(!std::isnan(intRad.real()) && !std::isnan(intRad.imag()));
    return leng * srcRWG->leng * intRad;
}
