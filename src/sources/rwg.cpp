#include "rwg.h"

RWG::RWG(
    std::shared_ptr<Excitation::PlaneWave> Einc,
    size_t rwgIdx,
    const Eigen::Vector4i& idxs,
    const TriVec& triangles)
    : Source(std::move(Einc), rwgIdx),
      tris({triangles[idxs[2]], triangles[idxs[3]]}),
      idx0(idxs[0]), idx1(idxs[1]),
      X0(Triangle::glVerts[idx0]), X1(Triangle::glVerts[idx1]),
      center((X0+X1)/2.0), 
      leng((X0-X1).norm()),
      bc(std::make_unique<BC>(this))
{
    // Find non-common vertices
    for (int i = 0; i < 2; ++i)
        for (const auto& iVerts : tris[i]->iVerts)
            if (iVerts != idxs[0] && iVerts != idxs[1]) {
                // idxpm[i] = iVerts;
                Xpm[i] = Triangle::glVerts[iVerts];
            }

    buildSubRWGs();

    buildVoltage(); // needs Xpm initialized!

    /*std::cout << '(' << X0 << ") (" << X1 << ") ("
        << Xpm[0] << ") (" << Xpm[1] << ") " << leng << '\n';*/
};

RWG::RWG(
    std::shared_ptr<Triangle> tri0,
    std::shared_ptr<Triangle> tri1)
    : tris( {std::move(tri0), std::move(tri1)} ) 
{ 
    // Find common and non-common vertices


}

void RWG::buildSubRWGs() {
    const Eigen::MatrixXi subIdxs{
        {0, 0, 4, 1, 2, 3},
        {5, 1, 5, 2, 3, 4}
    };

    int iTri = 0;
    TriVec centerSubtris;

    for (const auto& tri : tris) {
        const auto& subtris = tri->getSubtris({ Xpm[iTri], X0, X1 });
        //const auto& subtris = (!iTri ?
        //    tri->getSubtris({Xpm[iTri], X0, X1}) :
        //    tri->getSubtris({Xpm[iTri], X1, X0})
        //    );
        const int iTriBy8 = 8*iTri;
        
        for (int iSub = 0; iSub < 6; ++iSub) {
            auto subrwg = 
                std::make_shared<RWG>(subtris[subIdxs(0,iSub)], subtris[subIdxs(1,iSub)]);

            subrwgs[iTriBy8+iSub] = std::move(subrwg);
                
            if (iSub == 2 || iSub == 3)
                centerSubtris.push_back(subtris[iSub]);
        }
        /*subrwgs[iShift] = std::make_shared<RWG>(subtris[0], subtris[5]);
        subrwgs[iShift+1] = std::make_shared<RWG>(subtris[0], subtris[1]);
        subrwgs[iShift+2] = std::make_shared<RWG>(subtris[4], subtris[5]);
        subrwgs[iShift+3] = std::make_shared<RWG>(subtris[1], subtris[2]);
        subrwgs[iShift+4] = std::make_shared<RWG>(subtris[2], subtris[3]);
        subrwgs[iShift+5] = std::make_shared<RWG>(subtris[3], subtris[4]);*/

        ++iTri;
    }
    
    subrwgs[6] = std::make_shared<RWG>(centerSubtris[0], centerSubtris[2]);
    subrwgs[7] = std::make_shared<RWG>(centerSubtris[1], centerSubtris[3]);

    //for (int iSub = 0; iSub < 14; ++iSub)
    //    std::cout << subrwgs[iSub]->X0 << ' ' << subrwgs[iSub]->X1 << '\n';
}

vec3cd RWG::getIntegratedPlaneWave(const vec3d& kvec, bool doNumeric) const {
    using namespace Math;

    vec3cd rad = vec3cd::Zero();
    int triIdx = 0;

    if (doNumeric) {
        for (const auto& tri : tris) {
            for (const auto& [node, weight] : tri->getQuads())
                rad += weight * exp(iu*kvec.dot(node))
                        * (node - Xpm[triIdx])
                        * sign(triIdx);
            ++triIdx;
        }

        return leng * rad;
    }

    for (const auto& tri : tris) {
        const auto& Xs = tri->Xs, Ds = tri->Ds;

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
                f2 = (expI_alpha*(alphasq + 2.0*iu*alpha - 2.0) + 2.0) / (2.0*alpha*alphasq);

            radVec = -f1_alpha * (Xs[0] - Xpm[triIdx]) - iu*f2 * (Ds[0] - Ds[2]);

        } else {
            const cmplx
                I0 = (f0_alpha - f0_beta) / gamma,
                I1 = iu * (I0 + f1_alpha),
                I2 = -iu * (I0 + f1_beta);

            radVec = I0 * (Xs[0] - Xpm[triIdx]) + (I1*Ds[0] - I2*Ds[2]) / gamma;
        }

        rad += exp(iu*kvec.dot(Xs[0])) * radVec * sign(triIdx++);

        //if (rad.norm() > 1.0E3)
        //    std::cout << alpha << ' ' << beta << ' ' << gamma << ' ' << rad.norm() << '\n';
    }

    return leng * rad;
}

/* getIntegratedRad(src)
 * Return the radiated field due to src tested with this RWG
 */
cmplx RWG::getIntegratedRad(const std::shared_ptr<Source> src) const {

    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    const double k = Einc->wavenum;

    if (center == srcRWG->center) return 0.0; // TODO: Handle self-interactions
 
    cmplx intRad = 0.0;

    int obsTriIdx = 0;
    for (const auto& obsTri : tris) {

        const auto& obsQuads = obsTri->quads;
        const auto& obsXpm = Xpm[obsTriIdx];

        int srcTriIdx = 0;
        for (const auto& srcTri : srcRWG->tris) {

            const auto& srcQuads = srcTri->quads;
            const auto& srcXpm = srcRWG->Xpm[srcTriIdx];
            
            if (obsTri == srcTri) {
                ++srcTriIdx;
                continue; // TODO: Handle coincident tris
            }

            for (const auto& [obs, obsWeight] : obsQuads) {
                for (const auto& [src, srcWeight] : srcQuads) {

                    const vec3cd& rad = 
                        srcWeight 
                        * Math::dyadicG(obs-src, k) * (src-srcXpm) 
                        * Math::sign(srcTriIdx);

                    intRad += 
                        obsWeight 
                        * conj(rad.dot(obs-obsXpm)) // Hermitian dot!
                        * Math::sign(obsTriIdx);
                }
            }

            ++srcTriIdx;
        }

        ++obsTriIdx;
    }

    return leng * srcRWG->leng * intRad;
}