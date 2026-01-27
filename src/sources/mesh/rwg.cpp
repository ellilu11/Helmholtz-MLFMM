#include "rwg.h"

Mesh::RWG::RWG(
    std::shared_ptr<Excitation::PlaneWave> Einc, size_t iSrc, const vec4i& idx4)
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
}

void Mesh::RWG::refineRWGs() {
    triToSubs.resize(glTris.size());

    int iSub = 0;
    for (const auto& [edge, iTris] : fineEdgeToTri) {
        const vec4i& idx4 =
            { edge.first, edge.second, iTris[0], iTris[1] };

        if (iTris[1] < 0) continue; // Ignore edges not shared by two tris
        //std::cout << "Building SubRWG for edge ("
        //    << edge.first << ',' << edge.second << ") with tris # "
        //    << iTri[0] << ' ' << iTri[1] << '\n';

        glSubrwgs.emplace_back(iSub, idx4);
        fineEdgeToSub.emplace(edge, iSub);
        ++iSub;
    }

    //for (const auto& [edge, iSub] : fineEdgeToSub) {
    //    std::cout << "Edge (" << edge.first << ',' << edge.second
    //        << ") maps to subRWG # " << iSub << '\n';
    //}
}

vec3cd Mesh::RWG::getIntegratedPlaneWave(const vec3d& kvec, bool doNumeric) const {
    using namespace Math;

    const auto& tris = getTris();
    const auto& Xnc = getVertsNC();

    vec3cd rad = vec3cd::Zero();
    int iTri = 0;

    if (doNumeric) {
        for (const auto& tri : tris) {
            for (const auto& [node, weight] : tri.getQuads())
                rad += weight * exp(iu*kvec.dot(node))
                        * (node - Xnc[iTri])
                        * sign(iTri);
            ++iTri;
        }

        return leng * rad;
    }

    for (const auto& tri : tris) {
        const auto& Xs = tri.getVerts(), Ds = tri.Ds; // TODO: Compute Ds[0] and Ds[2] here

        const double alpha = kvec.dot(Ds[0]), beta = -kvec.dot(Ds[2]), gamma = alpha-beta;

        const double alphasq = alpha*alpha, betasq = beta*beta;

        const cmplx expI_alpha = exp(iu*alpha), expI_beta = exp(iu*beta);

        const cmplx // TODO: Only compute if gamma != 0
            f0_alpha = (approxZero(alpha) ? -iu : (1.0 - expI_alpha) / alpha),
            f0_beta = (approxZero(beta) ? -iu : (1.0 - expI_beta) / beta);

        const cmplx
            f1_alpha = (approxZero(alpha) ? -0.5 : (1.0 - (1.0 - iu*alpha) * expI_alpha) / alphasq),
            f1_beta = (approxZero(beta) ? -0.5 : (1.0 - (1.0 - iu*beta) * expI_beta) / betasq);

        vec3cd radVec;

        if (approxZero(gamma)) {
            // TODO: Handle alpha -> 0
            const cmplx
                f2 = (approxZero(alpha) ? 
                        iu/6.0 : 
                        (expI_alpha*(alphasq + 2.0*iu*alpha - 2.0) + 2.0) / (2.0*alpha*alphasq)
                     );

            radVec = -f1_alpha * (Xs[0] - Xnc[iTri]) - iu*f2 * (Ds[0] - Ds[2]);

        } else {
            const cmplx
                I0 = (f0_alpha - f0_beta) / gamma,
                I1 = iu * (I0 + f1_alpha),
                I2 = -iu * (I0 + f1_beta);

            radVec = I0 * (Xs[0] - Xnc[iTri]) + (I1*Ds[0] - I2*Ds[2]) / gamma;
        }

        rad += exp(iu*kvec.dot(Xs[0])) * radVec * sign(iTri++);
    }

    assert(!std::isnan(rad.norm()));
    return leng * rad;
}

/* getIntegratedRad(src)
 * Return the radiated field due to src tested with this RWG
 */
cmplx Mesh::RWG::getIntegratedRad(const std::shared_ptr<Source> src) const {

    const auto srcRWG = dynamic_pointer_cast<RWG>(src);

    const auto& tris = getTris();
    const auto& vertsNC = getVertsNC();

    const auto& srcTris = srcRWG->getTris();
    const auto& srcVertsNC = srcRWG->getVertsNC();

    const double k = Einc->wavenum;

    if (iSrc == src->getIdx()) return 0.0; // TODO: Handle self-interactions
 
    cmplx intRad = 0.0;

    int obsPairIdx = 0;
    for (const auto& obsTri : tris) {

        const auto& obsQuads = obsTri.quads;
        const auto& obsVert = vertsNC[obsPairIdx];

        int srcPairIdx = 0;
        for (const auto& srcTri : srcTris) {

            const auto& srcQuads = srcTri.quads;
            const auto& srcVert = srcVertsNC[srcPairIdx];
            
            if (obsTri.iTri == srcTri.iTri) {
                ++srcPairIdx;
                continue; // TODO: Handle coincident tris
            }

            for (const auto& [obs, obsWeight] : obsQuads) {
                for (const auto& [src, srcWeight] : srcQuads) {

                    const vec3cd& rad = 
                        srcWeight 
                        * Math::dyadicG(obs-src, k) * (src-srcVert) 
                        * Math::sign(srcPairIdx);

                    intRad += 
                        obsWeight 
                        * conj(rad.dot(obs-obsVert)) // Hermitian dot!
                        * Math::sign(obsPairIdx);
                }
            }

            ++srcPairIdx;
        }

        ++obsPairIdx;
    }

    assert(!std::isnan(intRad.real()) && !std::isnan(intRad.imag()));
    return leng * srcRWG->leng * intRad;
}