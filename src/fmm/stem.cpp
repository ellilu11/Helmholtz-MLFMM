#include "stem.h"

FMM::Stem::Stem(
    const SrcVec& srcs,
    const int branchIdx,
    Stem* const base)
    : Node(srcs, branchIdx, base)
{
    // Assign every src to a branch based on src center relative to node center
    std::array<SrcVec, 8> branchSrcs;

    for (const auto& src : srcs)
        branchSrcs[Math::bools2Idx(src->getCenter() > center)].push_back(src);

    // Construct branch nodes
    for (size_t k = 0; k < branchSrcs.size(); ++k) {
        std::shared_ptr<Node> branch;

        if (branchSrcs[k].size() > config.maxNodeSrcs)
            branch = std::make_shared<Stem>(branchSrcs[k], k, this);
        else
            branch = std::make_shared<Leaf>(branchSrcs[k], k, this);

        branches.push_back(std::move(branch));
    }
}

/* buildNeighbors()
 * Find all neighbor nodes of equal or greater size
 */
void FMM::Stem::buildNeighbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);

        if (nbor) nbors.push_back(nbor);
    }

    assert(nbors.size() <= numDir);
}

/* buildLists()
 * Find neighbor and interaction lists.
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void FMM::Stem::buildLists() {
    if (!isRoot()) {
        buildNeighbors();
        buildInteractionList();
        pushSelfToNearNonNbors();
    }

    for (const auto& branch : branches)
        branch->buildLists();
}

void FMM::Stem::resizeCoeffs() {
    const auto [nth, nph] = angles[level].getNumAngles();

    coeffs.resize(nth*nph);
    localCoeffs.resize(nth*nph);

    for (const auto& branch : branches)
        branch->resizeCoeffs();
}

/* buildMpoleCoeffs()
 * (M2M) Build mpole coeffs by merging branch mpole coeffs
 */
void FMM::Stem::buildMpoleCoeffs() {
    const int order = config.interpOrder;

    const auto [mth, mph] = angles[level+1].getNumAngles();

    coeffs.fillZero();

    for (const auto& branch : branches) {
        if (branch->isSrcless()) continue;

        branch->buildMpoleCoeffs();
        const auto& branchCoeffs = branch->getMpoleCoeffs();

        auto start = Clock::now();

        // Shift branch coeffs to center of this node
        const auto& dX = center - branch->getCenter();

        Coeffs shiftedCoeffs(mth*mph);
        for (int iDir = 0; iDir < mth*mph; ++iDir) {
            const vec3d& kvec = angles[level+1].khat[iDir] * wavenum;
            cmplx shift = exp(iu*kvec.dot(dX));

            shiftedCoeffs.theta[iDir] = shift * branchCoeffs.theta[iDir];
            shiftedCoeffs.phi[iDir] = shift * branchCoeffs.phi[iDir];
        }

        // Interpolate shifted coeffs to this node's angular grid
        addInterpCoeffs(shiftedCoeffs, coeffs, level+1, level);

        t.M2M += Clock::now() - start;
    }
}

/* getShiftedLocalCoeffs(branchIdx)
 * (L2L) Return local coeffs shifted to center of branch labeled by branchIdx
 * branchIdx : index of branch \in {0, ..., 7}
 */
FMM::Coeffs FMM::Stem::getShiftedLocalCoeffs(int branchIdx) const {

    const auto [mth, mph] = angles[level].getNumAngles();
    const auto [nth, nph] = angles[level+1].getNumAngles();

    Coeffs outCoeffs(nth*nph, 0.0);
    if (iList.empty()) return outCoeffs;

    // Shift local coeffs to center of branch
    const auto& dX = branches[branchIdx]->getCenter() - center;

    Coeffs shiftedCoeffs(mth*mph);
    for (int iDir = 0; iDir < mth*mph; ++iDir) {
        const vec3d& kvec = angles[level].khat[iDir] * wavenum;
        cmplx shift = exp(iu*kvec.dot(dX));

        shiftedCoeffs.theta[iDir] = shift * localCoeffs.theta[iDir];
        shiftedCoeffs.phi[iDir] = shift * localCoeffs.phi[iDir];
    }

    // Anterpolate shifted coeffs to branch's angular grid
    addAnterpCoeffs(shiftedCoeffs, outCoeffs, level, level+1);

    return outCoeffs;
}

void FMM::Stem::addInterpCoeffs(
    const Coeffs& inCoeffs, Coeffs& outCoeffs, int srcLvl, int tgtLvl)
{
    const int order = config.interpOrder;

    const auto [mth, mph] = angles[srcLvl].getNumAngles();
    const auto [nth, nph] = angles[tgtLvl].getNumAngles();
    assert(!(mph%2)); // mph needs to be even

    const int tblLvl = std::min(srcLvl, tgtLvl);
    const auto& interpTheta = tables[tblLvl].interpTheta;
    const auto& interpPhi = tables[tblLvl].interpPhi;

    // Interpolate over theta
    Coeffs innerCoeffs(nth*mph);

    for (int jth = 0; jth < nth; ++jth) {
        const auto [interp, nearIdx] = interpTheta[jth];

        for (int ith = nearIdx+1-order, k = 0; ith <= nearIdx+order; ++ith, ++k) {
            const int ith_flipped = Math::flipIdxToRange(ith, mth);

            const bool outOfRange = ith != ith_flipped; // jth < 0 || jth >= mth;

            for (int iph = 0; iph < mph; ++iph) {
                int iph_shifted = iph;
                if (outOfRange) iph_shifted += ((iph < mph/2) ? mph/2 : -mph/2);

                size_t idxInner = jth*mph+iph, idxIn = ith_flipped*mph+iph_shifted;

                innerCoeffs.theta[idxInner] += 
                    interp[k] * inCoeffs.theta[idxIn] * Math::sign(outOfRange); 
                innerCoeffs.phi[idxInner] += 
                    interp[k] * inCoeffs.phi[idxIn] * Math::sign(outOfRange);
            }
        }
    }

    // Interpolate over phi
    for (int jph = 0; jph < nph; ++jph) {
        const auto [interp, nearIdx] = interpPhi[jph];

        for (int iph = nearIdx+1-order, k = 0; iph <= nearIdx+order; ++iph, ++k) {
            const int iph_wrapped = Math::wrapIdxToRange(iph, mph);

            for (int jth = 0; jth < nth; ++jth) {
                size_t idxOut = jth*nph+jph, idxInner = jth*mph+iph_wrapped;

                outCoeffs.theta[idxOut] += interp[k] * innerCoeffs.theta[idxInner];
                outCoeffs.phi[idxOut] += interp[k] * innerCoeffs.phi[idxInner];
            }
        }
    }
}

void FMM::Stem::addAnterpCoeffs(
    const Coeffs& inCoeffs, Coeffs& outCoeffs, int srcLvl, int tgtLvl)
{
    const int order = config.interpOrder;

    const auto [mth, mph] = angles[srcLvl].getNumAngles();
    const auto [nth, nph] = angles[tgtLvl].getNumAngles();
    assert(!(nph%2)); // nph needs to be even

    const int tblLvl = std::min(srcLvl, tgtLvl);
    const auto& interpTheta = tables[tblLvl].interpTheta;
    const auto& interpPhi = tables[tblLvl].interpPhi;

    // Anterpolate over extended phi
    Coeffs innerCoeffs(mth*nph);
    for (int iph = 0; iph < mph; ++iph) {
        const auto [interp, nearIdx] = interpPhi[iph];

        for (int jph = -order; jph < nph+order; ++jph) {
            const int k = jph - (nearIdx+1-order);

            // If iph \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            const int jph_wrapped = Math::wrapIdxToRange(jph,nph);

            for (int ith = 0; ith < mth; ++ith) {
                size_t idxInner = ith*nph+jph_wrapped, idxIn = ith*mph+iph;

                innerCoeffs.theta[idxInner] += interp[k] * inCoeffs.theta[idxIn];
                innerCoeffs.phi[idxInner] += interp[k] * inCoeffs.phi[idxIn];
            }
        }
    }

    // Anterpolate over extended theta
    for (int ith = 0; ith < mth; ++ith) {
        const auto [interp, nearIdx] = interpTheta[ith];

        for (int jth = -order; jth < nth+order; ++jth) {  
            const int k = jth - (nearIdx+1-order);

            // If ith \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            const int jth_flipped = Math::flipIdxToRange(jth, nth);

            const bool outOfRange = jth != jth_flipped; // jth < 0 || jth >= nth;

            for (int jph = 0; jph < nph; ++jph) {
                int jph_shifted = jph;
                if (outOfRange) jph_shifted += ((jph < nph/2) ? nph/2 : -nph/2);

                size_t idxOut = jth_flipped*nph+jph_shifted, idxInner = ith*nph+jph;

                outCoeffs.theta[idxOut] +=
                    interp[k] * innerCoeffs.theta[idxInner] * Math::sign(outOfRange);
                outCoeffs.phi[idxOut] +=
                    interp[k] * innerCoeffs.phi[idxInner] * Math::sign(outOfRange);
            }
        }
    }
}

/* buildLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void FMM::Stem::buildLocalCoeffs() {
    if (!isRoot()) {
        auto start = Clock::now();
        translateCoeffs();
        t.M2L += Clock::now() - start;

        evalLeafIlistSols();

        start = Clock::now();
        if (!base->isRoot()) {
            localCoeffs = localCoeffs
                + dynamic_cast<Stem*>(base)->getShiftedLocalCoeffs(branchIdx);
        }

        t.L2L += Clock::now() - start;
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}
