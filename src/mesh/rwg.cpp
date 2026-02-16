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
        const auto& Xs = tri.getVerts(), Ds = tri.Ds;
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
            const cmplx
                f2 = (approxZero(alpha) ? iu/6.0 : 
                      (expI_alpha*(alphasq + 2.0*iu*alpha - 2.0) + 2.0) / (2.0*alpha*alphasq) );
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
    double k = Einc->wavenum;
 
    cmplx intRad = 0.0, intRadFull = 0.0;

    int iObsTri = 0;
    for (const auto& obsTri : getTris()) {
        const auto& obsNC = getVertsNC()[iObsTri];

        int iSrcTri = 0;
        for (const auto& srcTri : srcRWG->getTris()) {
            const auto& srcNC = srcRWG->getVertsNC()[iSrcTri];
            const auto& srcNCproj = srcTri.projectToPlane(srcNC);

            const int nCommon = obsTri.getNumCommonVerts(srcTri);

            for (const auto& [obs, obsWeight] : obsTri.triQuads) {
                if (!nCommon) { // No common vertices
                    for (const auto& [src, srcWeight] : srcTri.triQuads) {
                        const double r = (obs-src).norm();
                        const cmplx G = exp(iu*k*r) / r;
                        intRad += ((obs-obsNC).dot(src-srcNC) - 4.0 / (k*k)) * G
                            * obsWeight * srcWeight
                            * Math::sign(iObsTri) * Math::sign(iSrcTri);

                        //
                        const auto& dyadic = Math::dyadicG(obs-src, k);
                        intRadFull +=
                            (obs-obsNC).dot(dyadic*(src-srcNC))
                            * obsWeight * srcWeight
                            * Math::sign(iObsTri) * Math::sign(iSrcTri);
                        //
                    }
                } else { // Common vertex or common edge
                    // Add contribution from (e^(ikR)-1)/R term (numerically)
                    for (const auto& [src, srcWeight] : srcTri.triQuads) {
                        const double r = (obs-src).norm();
                        const cmplx G = (exp(iu*k*r)-1.0) / r;
                        intRad += ((obs-obsNC).dot(src-srcNC) - 4.0 / (k*k)) * G
                            * obsWeight * srcWeight
                            * Math::sign(iObsTri) * Math::sign(iSrcTri);
                    }

                    // Add contribution from 1/R term (analytically)
                    const auto& [scaRad, vecRad, obsProj] = srcTri.getNearIntegrated(obs);
                    intRad +=
                        ((obs-obsNC).dot(vecRad+(obsProj-srcNCproj)*scaRad)
                            - 4.0/(k*k)*scaRad)
                        * obsWeight * Math::sign(iObsTri) * Math::sign(iSrcTri);
                }
            }

            // For common triangles, add contribution from 1/R term (analytically)
            if (obsTri.iTri == srcTri.iTri) {
                const auto [V0, V1, V2] = obsTri.getVerts();
                double a00 = V0.dot(V0), a01 = V0.dot(V1), a02 = V0.dot(V2); // cache?

                intRad += obsTri.selfInts[0]
                    + obsTri.selfInts[1] * (-2.0*a00 + 2.0*a01 + (V0-V1).dot(srcNC+obsNC))
                    + obsTri.selfInts[2] * (-2.0*a00 + 2.0*a02 + (V0-V2).dot(srcNC+obsNC))
                    + obsTri.selfInts[3] * (a00 - V0.dot(srcNC+obsNC) + srcNC.dot(obsNC) - 4.0/(k*k));
            }

            /* Using precomputed moments
            const auto [mm0, mm1, mm2, mm3] =
                Mesh::glRadMoments.at(makeUnordered(obsTri.iTri, srcTri.iTri));
            intRad +=
                (mm0 - obsNC.dot(mm1) - mm2.dot(srcNC) + mm3 * (obsNC.dot(srcNC) - 4.0/(k*k)))
                * Math::sign(iObsTri) * Math::sign(iSrcTri);
            */

            ++iSrcTri;
        }

        ++iObsTri;
    }

    std::cout << std::setprecision(15)
        << "intRad (bilinear): " << intRad << '\n'
        << "intRadFull (full): " << intRadFull << '\n';

    //std::cout << std::setprecision(15) 
    //    << (intRad.real()-intRadFull.real()) / intRadFull.real() << ' '
    //    << (intRad.imag()-intRadFull.imag()) / intRadFull.imag() << '\n';

    assert(!std::isnan(intRad.real()) && !std::isnan(intRad.imag()));
    return leng * srcRWG->leng * intRad;
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
