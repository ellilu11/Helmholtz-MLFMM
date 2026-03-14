#include "farfield.h"

FMM::Farfield::Farfield(const std::shared_ptr<FMM::Node>& root) {
    if (root->isLeaf()) return;

    buildLevels();

    buildRadPats();

    resizeCoeffs(root);
}

/* buildLevels()
 * Build angular samples and interpolation/translation tables for each level
 */
void FMM::Farfield::buildLevels() {
    std::cout << " Building FMM operators...        ";

    auto start = Clock::now();

    angles.reserve(maxLevel+1);
    for (int level = 0; level <= maxLevel; ++level)
        angles.emplace_back(level);

    Tables::buildDists();
    tables.reserve(maxLevel+1);
    for (int level = 0; level <= maxLevel; ++level)
        tables.emplace_back(level, maxLevel);
    Tables::clearDists();

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";
}

void FMM::Farfield::buildRadPats() {
    std::cout << " Building plane wave expansions...";

    auto start = Clock::now();

    for (const auto& leaf : leaves) 
        buildRadPatsOfNode(leaf);

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";
}

/* resizeCoeffs(node)
 * Resize mpole and local coeffs of node and its branches to match 
 * number of directions at each level
 */
void FMM::Farfield::resizeCoeffs(const std::shared_ptr<FMM::Node>& node) {
    size_t nDir = angles[node->level].getNumDirs();
    node->coeffs.resize(nDir);
    node->localCoeffs.resize(nDir);

    for (const auto& branch : node->branches)
        resizeCoeffs(branch);
}

/* buildRadPats()
 * Build radiation patterns due to sources in leaf node
 */
void FMM::Farfield::buildRadPatsOfNode(const std::shared_ptr<FMM::Node>& node) {
    if (node->isSrcless()) return;
    const Angles& angles_lvl = angles[node->level];
    size_t nDir = angles_lvl.getNumDirs();

    node->radPats.resize(node->srcs.size());
    size_t iSrc = 0;
    for (const auto& src : node->srcs) {
        Coeffs radPat(nDir);

        for (int iDir = 0; iDir < nDir; ++iDir) {
            vec3d kvec = angles_lvl.khat[iDir] * config.k;
            const mat23d& toThPh = angles_lvl.toThPh[iDir];

            radPat.setCoeffAlongDir(
                toThPh * src->getRadAlongDir(node->center, kvec), iDir);
        }

        node->radPats[iSrc++] = std::move(radPat);
    }
}

/* buildMpoleCoeffs()
 * (S2M) Build multipole coefficients from sources in this leaf
 */
void FMM::Farfield::buildMpoleCoeffs(const std::shared_ptr<FMM::Node>& node) {
    Coeffs& coeffs = node->coeffs;

    coeffs.fillZero();
    if (node->isSrcless() || node->isRoot()) return;

    auto start = Clock::now();

    size_t nDir = angles[node->level].getNumDirs();
    Eigen::Map<arrXcd> coeffsTheta(coeffs.theta.data(), nDir);
    Eigen::Map<arrXcd> coeffsPhi(coeffs.phi.data(), nDir);

    size_t iSrc = 0;
    for (const auto& src : node->srcs) {
        Coeffs& radPat = node->radPats[iSrc++];

        Eigen::Map<arrXcd> radPatTheta(radPat.theta.data(), nDir);
        Eigen::Map<arrXcd> radPatPhi(radPat.phi.data(), nDir);

        coeffsTheta += Solver::lvec[src->getIdx()] * radPatTheta;
        coeffsPhi += Solver::lvec[src->getIdx()] * radPatPhi;
    }

    t.S2M += Clock::now() - start;
}

/* mergeMpoleCoeffs()
 * (M2M) Build mpole coeffs by merging branch mpole coeffs
 */
void FMM::Farfield::mergeMpoleCoeffs(const std::shared_ptr<FMM::Node>& node) {
    int order = config.interpOrder;
    int level = node->level;
    size_t mDir = angles[level+1].getNumDirs();

    node->coeffs.fillZero();

    for (const auto& branch : node->branches) {
        if (branch->isSrcless()) continue;

        if (branch->isLeaf()) buildMpoleCoeffs(branch);
        else mergeMpoleCoeffs(branch);
        const Coeffs& branchCoeffs = branch->coeffs;

        auto start = Clock::now();

        // Shift branch coeffs to center of this node
        vec3d dX = node->center - branch->center;

        Coeffs shiftedCoeffs(mDir);
        for (int iDir = 0; iDir < mDir; ++iDir) {
            const vec3d& kvec = angles[level+1].khat[iDir] * config.k;
            cmplx shift = exp(iu*kvec.dot(dX));

            shiftedCoeffs.theta[iDir] = shift * branchCoeffs.theta[iDir];
            shiftedCoeffs.phi[iDir] = shift * branchCoeffs.phi[iDir];
        }

        // Interpolate shifted coeffs to this node's angular grid
        addInterpCoeffs(shiftedCoeffs, node->coeffs, level+1, level);

        t.M2M += Clock::now() - start;
    }
}

/* translateCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center,
 * then apply integration weights for anterpolation
 */
void FMM::Farfield::translateCoeffs(const std::shared_ptr<FMM::Node>& node) {
    if (node->iList.empty()) return;

    Coeffs& localCoeffs = node->localCoeffs;
    localCoeffs.fillZero();
    size_t nDir = localCoeffs.size();

    // Translate mpole coeffs into local coeffs
    const VecHashMap<arrXcd>& transl = tables[node->level].transl;
    Eigen::Map<arrXcd> localTheta(localCoeffs.theta.data(), nDir);
    Eigen::Map<arrXcd> localPhi(localCoeffs.phi.data(), nDir);

    for (const auto& srcNode : node->iList) {
        vec3d dX = node->center - srcNode->center;
        arrXcd transl_dX = transl.at(dX/node->nodeLeng);

        Eigen::Map<arrXcd> mpoleTheta(srcNode->coeffs.theta.data(), nDir);
        Eigen::Map<arrXcd> mpolePhi(srcNode->coeffs.phi.data(), nDir);

        localTheta += transl_dX * mpoleTheta;
        localPhi += transl_dX * mpolePhi;
    }

    // Apply integration weights
    const Angles& angles_lvl = angles[node->level];
    auto [nth, nph] = angles_lvl.getNumAngles();
    double phiWeight = 2.0*PI / static_cast<double>(nph);

    size_t iDir = 0;
    for (int ith = 0; ith < nth; ++ith) {
        double thetaWeight = angles_lvl.weights[ith];

        for (int iph = 0; iph < nph; ++iph) {
            localCoeffs.theta[iDir] *= thetaWeight * phiWeight;
            localCoeffs.phi[iDir] *= thetaWeight * phiWeight;
            ++iDir;
        }
    }
}

/* getShiftedLocalCoeffs(branchIdx)
 * (L2L) Return local coeffs shifted to center of branch labeled by branchIdx
 * branchIdx : index of branch \in {0, ..., 7}
 */
FMM::Coeffs FMM::Farfield::getShiftedLocalCoeffs(
    FMM::Node* node, int branchIdx) const 
{
    int level = node->level;
    size_t mDir = angles[level].getNumDirs();
    size_t nDir = angles[level+1].getNumDirs();

    Coeffs outCoeffs(nDir);
    if (node->iList.empty()) return outCoeffs;

    // Shift local coeffs to center of branch
    vec3d dX = node->branches[branchIdx]->center - node->center;

    Coeffs shiftedCoeffs(mDir);
    for (int iDir = 0; iDir < mDir; ++iDir) {
        vec3d kvec = angles[level].khat[iDir] * config.k;
        cmplx shift = exp(iu*kvec.dot(dX));

        shiftedCoeffs.theta[iDir] = shift * node->localCoeffs.theta[iDir];
        shiftedCoeffs.phi[iDir] = shift * node->localCoeffs.phi[iDir];
    }

    // Anterpolate shifted coeffs to branch's angular grid
    addAnterpCoeffs(shiftedCoeffs, outCoeffs, level, level+1);

    return outCoeffs;
}

/* buildLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void FMM::Farfield::buildLocalCoeffs(const std::shared_ptr<FMM::Node>& node) {
    if (!node->isRoot()) {
        // M2L
        auto start = Clock::now();
        translateCoeffs(node);
        t.M2L += Clock::now() - start;

        // Add L2L to M2L
        start = Clock::now();
        if (!node->base->isRoot()) // TODO: Avoid adding Coeffs directly
            node->localCoeffs += 
                getShiftedLocalCoeffs(node->base, node->branchIdx);
        t.L2L += Clock::now() - start;
    }

    for (const auto& branch : node->branches)
        buildLocalCoeffs(branch);
}

/* evalFarSols()
 * (L2T) Evaluate sols from local expansion due to far nodes
 */
void FMM::Farfield::evalFarSols(const std::shared_ptr<FMM::Node>& node) {
    int level = node->level;
    if (node->isSrcless() || level <= 1) return;

    size_t nDir = angles[level].getNumDirs();
    Eigen::Map<arrXcd> localTheta(node->localCoeffs.theta.data(), nDir);
    Eigen::Map<arrXcd> localPhi(node->localCoeffs.phi.data(), nDir);

    size_t iObs = 0;
    for (const auto& obs : node->srcs) {
        cmplx intRad = 0;

        Coeffs& radPat = node->radPats[iObs++];
        Eigen::Map<arrXcd> radPatTheta(radPat.theta.data(), nDir);
        Eigen::Map<arrXcd> radPatPhi(radPat.phi.data(), nDir);

        intRad += (radPatTheta.conjugate() * localTheta +
            radPatPhi.conjugate() * localPhi).sum();

        Solver::rvec[obs->getIdx()] += Phys::C * config.k * intRad;
    }
}

/* evaluateSols()
 * (L2T) Evaluate sols from local expansion due to far nodes for all leaves
 */
void FMM::Farfield::evaluateSols() {
    auto start = Clock::now();

    for (const auto& leaf : leaves) 
        evalFarSols(leaf);

    t.L2T += Clock::now() - start;
}