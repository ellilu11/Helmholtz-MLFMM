#include "leaf.h"

FMM::Leaf::Leaf(
    const SrcVec& srcs,
    const int branchIdx,
    Stem* const base)
    : Node(srcs, branchIdx, base), 
      leafPairIdx(0), nonNearPairIdx(0)
{
    maxLevel = std::max(level, maxLevel);

    /* Assign indices to all sources in this leaf
    for (const auto& src : srcs) {
        src->setIdx(glSrcIdx++);
        std::cout << src->getIdx() << ' ';
    }
    if (!isSrcless()) std::cout << '\n';
    */
}

/* buildNeighbors()
 * Find all neighbor nodes of equal or greater size
 * Also find all neighbor leaves of equal or lesser size (list 1)
 */
void FMM::Leaf::buildNeighbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);

        if (nbor) {
            nbors.push_back(nbor);
            auto nbors = getNeighborsLeqSize(nbor, dir);
            nearNbors.insert(nearNbors.end(), nbors.begin(), nbors.end());
        }
    }

    assert(nbors.size() <= numDir);
}

/* buildLists()
 * Add self to list of leaves 
 * Find neighbor and interaction lists
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void FMM::Leaf::buildLists() {
    leaves.push_back(shared_from_this()); 

    if (isRoot()) return;
    
    buildNeighbors();
    buildInteractionList();
    pushSelfToNearNonNbors();
}

void FMM::Leaf::resizeCoeffs() {
    const auto [nth, nph] = angles[level].getNumAngles();

    coeffs.resize(nth*nph);
    localCoeffs.resize(nth*nph);
}

/* findNearNborPairs()
 * From list of leaves, find all near neighbor leaf pairs
 */
void FMM::Leaf::findNearNborPairs() {

    for (const auto& leaf : leaves) {
        for (const auto& nbor : leaf->nearNbors) {
            auto nborLeaf = dynamic_pointer_cast<Leaf>(nbor);

            if (leaf < nborLeaf)
                nearPairs.emplace_back(leaf, nborLeaf);
        }
    }
}

void FMM::Leaf::buildNearRads() {
    findNearNborPairs();

    for (const auto& [obsLeaf, srcLeaf] : nearPairs) {
        assert(obsLeaf < srcLeaf);

        const size_t nObs = obsLeaf->srcs.size(), nSrcs = srcLeaf->srcs.size();

        auto leafPairRads = std::vector<cmplx>(nObs*nSrcs);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nObs; ++iObs) {
            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                const auto obs = obsLeaf->srcs[iObs], src = srcLeaf->srcs[iSrc];

                leafPairRads[pairIdx++] = obs->getIntegratedRad(src);
            }
        }

        obsLeaf->nearRads.push_back(leafPairRads);
    }

    for (const auto& [obsNode, srcNode] : nonNearPairs) {
        auto obsLeaf = dynamic_pointer_cast<Leaf>(obsNode);

        const size_t nObs = obsLeaf->srcs.size(), nSrcs = srcNode->getSrcs().size();

        auto nodePairRads = std::vector<cmplx>(nObs*nSrcs);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nObs; ++iObs) {
            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                const auto obs = obsLeaf->srcs[iObs], src = srcNode->getSrcs()[iSrc];

                nodePairRads[pairIdx++] = obs->getIntegratedRad(src);
            }
        }

        obsLeaf->nonNearRads.push_back(nodePairRads);
    }

    for (const auto& leaf : leaves) {
        for (size_t iObs = 0; iObs < leaf->srcs.size(); ++iObs) { // iObs = 0
            const auto& obs = leaf->srcs[iObs];

            for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
                const auto& src = leaf->srcs[iSrc];

                leaf->selfRads.push_back(obs->getIntegratedRad(src));
            }
        }
    }
}

/* buildRadPats()
 * Build radiation patterns due to sources in all leaves
 */
void FMM::Leaf::buildRadPats() {
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
FMM::Coeffs FMM::Leaf::buildMpoleCoeffs() {
    auto start = Clock::now();

    coeffs.fillZero();
    if (isSrcless() || isRoot()) return coeffs;

    size_t nDir = angles[level].getNumDirs();
    size_t iSrc = 0;

    for (const auto& src : srcs)
        coeffs += states.lvec[src->getIdx()] * radPats[iSrc++];
  
    t.S2M += Clock::now() - start;

    return coeffs;
}

/* buildLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void FMM::Leaf::buildLocalCoeffs() {
    if (isRoot()) return;

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

/* evalFarSols()
 * (L2T) Evaluate sols from local expansion due to far nodes
 */
void FMM::Leaf::evalFarSols() {
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

/* evalNearNonNborSols()
 * (M2T/S2T) Evaluate sols from mpole expansion due to list 3 nodes
 */
void FMM::Leaf::evalNearNonNborSols() {
    for (const auto& node : nearNonNbors)
        evalPairSols(node, nonNearRads[nonNearPairIdx++]);
    return;
}

/* (S2T) Evaluate sols at sources in this leaf due to sources in srcNode
 * and vice versa
 * srcNode : source node
 * rads    : precomputed radiation coefficients
 */
void FMM::Leaf::evalPairSols(
    const std::shared_ptr<Node> srcNode, const std::vector<cmplx>& rads) 
{
    const auto& srcSrcs = srcNode->getSrcs();

    const int nObs = srcs.size(), nSrcs = srcSrcs.size();

    std::vector<cmplx> solAtObss(nObs, 0.0);
    std::vector<cmplx> solAtSrcs(nSrcs, 0.0);

    int pairIdx = 0;
    for (size_t iObs = 0; iObs < nObs; ++iObs) {
        for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
            const auto obs = srcs[iObs], src = srcSrcs[iSrc];

            const cmplx rad = rads[pairIdx++];

            solAtObss[iObs] += states.lvec[src->getIdx()] * rad;
            solAtSrcs[iSrc] += states.lvec[obs->getIdx()] * rad;
        }
    }

    for (int n = 0; n < nObs; ++n)
        states.rvec[srcs[n]->getIdx()] += Phys::C * k * solAtObss[n];

    for (int n = 0; n < nSrcs; ++n)
        states.rvec[srcSrcs[n]->getIdx()] += Phys::C * k * solAtSrcs[n];
}

/* evalSelfSols()
 * (S2T) Evaluate sols at sources in this node due to other sources
 * in this node
 */
void FMM::Leaf::evalSelfSols() {
    const int nSrcs = srcs.size();

    std::vector<cmplx> solAtObss(nSrcs, 0.0);

    int pairIdx = 0;
    for (size_t iObs = 0; iObs < nSrcs; ++iObs) { // iObs = 0
        for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
            auto obs = srcs[iObs], src = srcs[iSrc];

            const cmplx rad = selfRads[pairIdx++];

            solAtObss[iObs] += states.lvec[src->getIdx()] * rad;
            solAtObss[iSrc] += states.lvec[obs->getIdx()] * rad;
        }
    }

    for (int n = 0; n < nSrcs; ++n)
        states.rvec[srcs[n]->getIdx()] += Phys::C * k * solAtObss[n];
}

/* evaluateSols()
 * Sum solutions at all sources in all leaves 
 */ 
void FMM::Leaf::evaluateSols() {
    auto start = Clock::now();
  
    for (const auto& [obsLeaf, srcLeaf] : nearPairs) {
        auto pairIdx = obsLeaf->leafPairIdx++;
        obsLeaf->evalPairSols(srcLeaf, obsLeaf->nearRads[pairIdx]);
    }

    for (const auto& leaf : leaves) {
        leaf->evalFarSols();

        leaf->evalNearNonNborSols();

        leaf->evalSelfSols();

        leaf->leafPairIdx = 0;
        leaf->nonNearPairIdx = 0;
    }

    t.L2T += Clock::now() - start;
}