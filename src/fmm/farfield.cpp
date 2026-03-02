#include "node.h"

/* buildRadPats()
 * Build radiation patterns due to sources in this leaf
 */
void FMM::Node::buildRadPats() {
    // assert(isLeaf());
    const Angles& angles_lvl = angles[level];
    size_t nDir = angles_lvl.getNumDirs();

    radPats.resize(srcs.size());
    size_t iSrc = 0;
    for (const auto& src : srcs) {
        Coeffs radPat(nDir);

        for (int iDir = 0; iDir < nDir; ++iDir) {
            vec3d kvec = angles_lvl.khat[iDir] * k;
            const mat23d& toThPh = angles_lvl.toThPh[iDir];

            radPat.setCoeffAlongDir(
                toThPh * src->getRadAlongDir(center, kvec), iDir);
        }

        radPats[iSrc++] = std::move(radPat);
    }
}

/* buildMpoleCoeffs()
 * (S2M) Build multipole coefficients from sources in this leaf
 */
FMM::Coeffs FMM::Node::buildMpoleCoeffs() {
    coeffs.fill(0.0);
    if (isSrcless() || isRoot()) return coeffs;

    auto start = Clock::now();

    size_t nDir = angles[level].getNumDirs();
    Eigen::Map<arrXcd> coeffsArr(coeffs.vals.data(), 2*nDir);

    size_t iSrc = 0;
    for (const auto& src : srcs) {
        Coeffs& radPat = radPats[iSrc++];
        Eigen::Map<arrXcd> radPatArr(radPat.vals.data(), 2*nDir);

        coeffsArr += states.lvec[src->getIdx()] * radPatArr;
    }

    t.S2M += Clock::now() - start;

    return coeffs;
}

/* mergeMpoleCoeffs()
 * (M2M) Build mpole coeffs by merging branch mpole coeffs
 */
FMM::Coeffs FMM::Node::mergeMpoleCoeffs() {
    int order = config.interpOrder;
    size_t mDir = angles[level+1].getNumDirs();

    coeffs.fill(0.0);

    for (const auto& branch : branches) {
        if (branch->isSrcless()) continue;

        const Coeffs& branchCoeffs = (branch->isLeaf() ? 
            branch->buildMpoleCoeffs() : branch->mergeMpoleCoeffs() );

        auto start = Clock::now();

        // Shift branch coeffs to center of this node
        vec3d dX = center - branch->getCenter();

        Coeffs shiftedCoeffs(mDir);
        for (int iDir = 0; iDir < mDir; ++iDir) {
            const vec3d& kvec = angles[level+1].khat[iDir] * k;
            cmplx shift = exp(iu*kvec.dot(dX));

            shiftedCoeffs.vals[iDir] = shift * branchCoeffs.vals[iDir];
            shiftedCoeffs.vals[mDir+iDir] = shift * branchCoeffs.vals[mDir+iDir];
        }

        // Interpolate shifted coeffs to this node's angular grid
        addInterpCoeffs(shiftedCoeffs, coeffs, level+1, level);

        t.M2M += Clock::now() - start;
    }

    return coeffs;
}

/* translateCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center,
 * then apply integration weights for anterpolation
 */
void FMM::Node::translateCoeffs() {
    if (iList.empty()) return;

    localCoeffs.fill(0.0);
    size_t nDir = localCoeffs.size();

    // Translate mpole coeffs into local coeffs
    const VecHashMap<arrXcd>& transl = tables[level].transl;
    Eigen::Map<arrXcd> localArr(localCoeffs.vals.data(), 2*nDir);

    for (const auto& node : iList) {
        vec3d dX = center - node->center;
        arrXcd transl_dX = transl.at(dX/nodeLeng);

        Eigen::Map<arrXcd> mpoleArr(node->coeffs.vals.data(), 2*nDir);

        localArr += transl_dX * mpoleArr;
    }

    // Apply integration weights
    const Angles& angles_lvl = angles[level];
    auto [nth, nph] = angles_lvl.getNumAngles();
    double phiWeight = 2.0*PI / static_cast<double>(nph);

    size_t iDir = 0;
    for (int ith = 0; ith < nth; ++ith) {
        double thetaWeight = angles_lvl.weights[ith];

        for (int iph = 0; iph < nph; ++iph) {
            localCoeffs.vals[iDir] *= thetaWeight * phiWeight;
            localCoeffs.vals[nDir+iDir] *= thetaWeight * phiWeight;
            ++iDir;
        }
    }
}

/* getShiftedLocalCoeffs(branchIdx)
 * (L2L) Return local coeffs shifted to center of branch labeled by branchIdx
 * branchIdx : index of branch \in {0, ..., 7}
 */
FMM::Coeffs FMM::Node::getShiftedLocalCoeffs(int branchIdx) const {
    size_t mDir = angles[level].getNumDirs();
    size_t nDir = angles[level+1].getNumDirs();

    Coeffs outCoeffs(nDir);
    if (iList.empty()) return outCoeffs;

    // Shift local coeffs to center of branch
    vec3d dX = branches[branchIdx]->getCenter() - center;

    Coeffs shiftedCoeffs(mDir);
    for (int iDir = 0; iDir < mDir; ++iDir) {
        vec3d kvec = angles[level].khat[iDir] * k;
        cmplx shift = exp(iu*kvec.dot(dX));

        shiftedCoeffs.vals[iDir] = shift * localCoeffs.vals[iDir];
        shiftedCoeffs.vals[mDir+iDir] = shift * localCoeffs.vals[mDir+iDir];
    }

    // Anterpolate shifted coeffs to branch's angular grid
    addAnterpCoeffs(shiftedCoeffs, outCoeffs, level, level+1);

    return outCoeffs;
}

/* buildLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void FMM::Node::buildLocalCoeffs() {
    if (!isRoot()) {
        // M2L
        auto start = Clock::now();
        translateCoeffs();
        t.M2L += Clock::now() - start;

        // Add L2L to M2L
        start = Clock::now();
        if (!base->isRoot())
            localCoeffs += base->getShiftedLocalCoeffs(branchIdx);
        t.L2L += Clock::now() - start;
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}

/* evalFarSols()
 * (L2T) Evaluate sols from local expansion due to far nodes
 */
void FMM::Node::evalFarSols() {
    if (isSrcless() || level <= 1) return;

    size_t nDir = angles[level].getNumDirs();
    Eigen::Map<arrXcd> localArr(localCoeffs.vals.data(), 2*nDir);

    size_t iObs = 0;
    for (const auto& obs : srcs) {
        cmplx intRad = 0;

        Coeffs& radPat = radPats[iObs++];
        Eigen::Map<arrXcd> radPatArr(radPat.vals.data(), 2*nDir);

        intRad += (radPatArr.conjugate() * localArr).sum();

        states.rvec[obs->getIdx()] += Phys::C * k * intRad;
    }
}