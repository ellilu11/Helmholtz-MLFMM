#include "node.h"

/* buildRadPats()
 * Build radiation patterns due to sources in all leaves
 */
void FMM::Node::buildRadPats() {
    for (const auto& leaf : leaves) {
        const auto& angles_lvl = angles[leaf->level];
        size_t nDir = angles_lvl.getNumDirs();

        leaf->radPats.resize(leaf->srcs.size());
        size_t iSrc = 0;
        for (const auto& src : leaf->srcs) {
            Coeffs radPat(nDir);

            for (int iDir = 0; iDir < nDir; ++iDir) {
                const auto& kvec = angles_lvl.khat[iDir] * k;
                const auto& toThPh = angles_lvl.toThPh[iDir];

                radPat.setCoeffAlongDir(
                    toThPh * src->getRadAlongDir(leaf->center, kvec), iDir);
            }

            leaf->radPats[iSrc++] = std::move(radPat);
        }
    }
}

/* buildMpoleCoeffs()
 * (S2M) Build multipole coefficients from sources in this node
 */
FMM::Coeffs FMM::Node::buildMpoleCoeffs() {
    auto start = Clock::now();

    coeffs.fillZero();
    if (isSrcless() || isRoot()) return coeffs;

    size_t iSrc = 0;
    for (const auto& src : srcs)
        coeffs += states.lvec[src->getIdx()] * radPats[iSrc++];

    t.S2M += Clock::now() - start;

    return coeffs;
}

/* mergeMpoleCoeffs()
 * (M2M) Build mpole coeffs by merging branch mpole coeffs
 */
FMM::Coeffs FMM::Node::mergeMpoleCoeffs() {
    int order = config.interpOrder;
    size_t mDir = angles[level+1].getNumDirs();

    coeffs.fillZero();

    for (const auto& branch : branches) {
        if (branch->isSrcless()) continue;

        const auto& branchCoeffs = (branch->isLeaf() ? 
            branch->buildMpoleCoeffs() : branch->mergeMpoleCoeffs() );

        auto start = Clock::now();

        // Shift branch coeffs to center of this node
        const auto& dX = center - branch->getCenter();

        Coeffs shiftedCoeffs(mDir);
        for (int iDir = 0; iDir < mDir; ++iDir) {
            const vec3d& kvec = angles[level+1].khat[iDir] * k;
            cmplx shift = exp(iu*kvec.dot(dX));

            shiftedCoeffs.theta[iDir] = shift * branchCoeffs.theta[iDir];
            shiftedCoeffs.phi[iDir] = shift * branchCoeffs.phi[iDir];
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

    localCoeffs.fillZero();
    size_t nDir = localCoeffs.size();

    // Translate mpole coeffs into local coeffs
    const VecHashMap<arrXcd>& transl = tables[level].transl;
    Eigen::Map<arrXcd> localTheta(localCoeffs.theta.data(), nDir);
    Eigen::Map<arrXcd> localPhi(localCoeffs.phi.data(), nDir);

    for (const auto& node : iList) {
        const auto& dX = center - node->center;
        const auto& transl_dX = transl.at(dX/nodeLeng);

        Eigen::Map<arrXcd> mpoleTheta(node->coeffs.theta.data(), nDir);
        Eigen::Map<arrXcd> mpolePhi(node->coeffs.phi.data(), nDir);

        localTheta += transl_dX * mpoleTheta;
        localPhi += transl_dX * mpolePhi;
    }

    // Apply integration weights
    const auto& angles_lvl = angles[level];
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
FMM::Coeffs FMM::Node::getShiftedLocalCoeffs(int branchIdx) const {
    size_t mDir = angles[level].getNumDirs();
    size_t nDir = angles[level+1].getNumDirs();

    Coeffs outCoeffs(nDir);
    if (iList.empty()) return outCoeffs;

    // Shift local coeffs to center of branch
    const auto& dX = branches[branchIdx]->getCenter() - center;

    Coeffs shiftedCoeffs(mDir);
    for (int iDir = 0; iDir < mDir; ++iDir) {
        const vec3d& kvec = angles[level].khat[iDir] * k;
        cmplx shift = exp(iu*kvec.dot(dX));

        shiftedCoeffs.theta[iDir] = shift * localCoeffs.theta[iDir];
        shiftedCoeffs.phi[iDir] = shift * localCoeffs.phi[iDir];
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

    Eigen::Map<arrXcd> localTheta(localCoeffs.theta.data(), nDir);
    Eigen::Map<arrXcd> localPhi(localCoeffs.phi.data(), nDir);

    size_t iObs = 0;
    for (const auto& obs : srcs) {
        cmplx intRad = 0;

        auto& radPat = radPats[iObs++];
        Eigen::Map<arrXcd> radPatTheta(radPat.theta.data(), nDir);
        Eigen::Map<arrXcd> radPatPhi(radPat.phi.data(), nDir);

        intRad += (radPatTheta.conjugate() * localTheta +
            radPatPhi.conjugate() * localPhi).sum();

        states.rvec[obs->getIdx()] += Phys::C * k * intRad;
    }
}

// TODO: Move into Farfield class
void FMM::Node::evaluateSols() {
    auto start = Clock::now();

    for (const auto& leaf : leaves)
        leaf->evalFarSols();

    t.L2T += Clock::now() - start;
}

void FMM::Node::printScatteredField(const std::string& fname, int nth, bool useMie) {
    namespace fs = std::filesystem;
    fs::path dir = "out/ff";
    std::error_code ec;

    std::cout << " Computing scattered farfield...\n";

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream farfile(dir/fname);
    farfile << std::setprecision(15) << std::scientific;

    // Also print out angles (coordinates of farsols)
    std::ofstream thfile(dir/"thetas.txt"); // phfile(dir/"phis.txt");
    thfile << std::setprecision(15) << std::scientific;

    double phi = 0.0; // pick phi = 0
    for (int ith = 0; ith < nth; ++ith) {
        double theta = (ith+0.5)*PI/static_cast<double>(nth);
        const auto& rhat = Math::fromSph(vec3d(1.0, theta, phi));

        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += states.currents[src->getIdx()] * src->getFarAlongDir(k*rhat);

        const vec2cd& far = Phys::C * k * Math::toThPh(theta, phi) * dirFar;

        farfile << norm(far[0]) << ' ' << norm(far[1]) << '\n'; // squared magnitude of theta and phi components
        thfile << theta << '\n';
    }
}
