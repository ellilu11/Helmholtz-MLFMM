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

/* getIntegratedPlaneWave(kvec, doNumeric)
 * Return integral of exp(ik dot r'} * f(r') dr' at this RWG
 */
vec3cd Mesh::RWG::getIntegratedPlaneWave(const vec3d& kvec, bool doNumeric) const {
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
        const auto& [scaRad, vecRad] = tri.getIntegratedPlaneWave(kvec);
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
        const vec3d& obsNC = getVertsNC()[iObsTri];

        int iSrcTri = 0;
        for (const auto& srcTri : srcRWG->getTris() ) {
            const vec3d& srcNC = srcRWG->getVertsNC()[iSrcTri];

            const auto& iTriPair = std::minmax(obsTri.iTri, srcTri.iTri);
            // std::cout << "(" << iTriPair.first << "," << iTriPair.second << ") ";

            const auto& triPair = glTriPairs.at(iTriPair);
            const auto& [m00, m10, m01, m11] = triPair.radMoments;

            cmplx pairRad = 
                m11 - m10.dot(srcNC) - obsNC.dot(m01) + (obsNC.dot(srcNC) - 4.0/(k*k))*m00;

            // For edge adjacent triangles, integrate 1/R term (analytically)
            if (triPair.nCommon >= 2)
                pairRad += triPair.getDoubleIntegratedInvR(obsNC, srcNC);

            /* For common triangles, add contribution from 1/R term (analytically)
            if (nCommon == 3)
                pairRad = obsTri.getDoubleSelfIntegratedInvR(obsNC, srcNC);
            */

            intRad += pairRad * Math::sign(iSrcTri) * Math::sign(iObsTri);
            ++iSrcTri;
        }

        ++iObsTri;
    }

    assert(!std::isnan(intRad.real()) && !std::isnan(intRad.imag()));
    return leng * srcRWG->leng * intRad;
}
