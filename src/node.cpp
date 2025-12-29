#include "node.h"

Config Node::config;
double Node::wavenum;
std::vector<realVec> Node::thetas;
std::vector<realVec> Node::thetaWeights;
std::vector<realVec> Node::phis;
std::vector<int> Node::Ls;
Tables Node::tables;

std::shared_ptr<vecXcd> Node::currents;
std::shared_ptr<vecXcd> Node::sols;

std::vector<NodePair> Node::nonNearPairs;

void Node::initNodes(
    const Config& config_, 
    const std::shared_ptr<Excitation::PlaneWave>& Einc,
    const std::unique_ptr<Solver>& solver) 
{
    config = config_;
    wavenum = Einc->wavenum;

    currents = std::move(solver->getQvec());
    sols = std::move(solver->getSols());
}

/* Node(particles,branchIdx,base)
 * particles : list of particles contained in this node
 * branchidx : index of this node relative to its base node
 * base      : pointer to base node
 */
Node::Node(
    const SrcVec& srcs,
    const int branchIdx,
    Node* const base)
    : srcs(srcs), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? config.rootLeng : base->nodeLeng / 2.0),
    level(base == nullptr ? 0 : base->level + 1),
    center(base == nullptr ? zeroVec :
        base->center + nodeLeng / 2.0 * Math::idx2pm(branchIdx)),
    nonNearPairIdx(0)
{
    numNodes++;
}

/* buildAngularSamples()
 * Compute theta and phi samples at each level
 */
void Node::buildAngularSamples() {

    constexpr double EPS_NR = 1.0E-9; // Newton-Raphson precision

    std::cout << "   (Lvl,Nth,Nph) =\n";

    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        const double nodeLeng = config.rootLeng / pow(2.0, lvl);

        // Use excess bandwidth formula
        const int tau = 
            ceil(
                (1.73*wavenum*nodeLeng +
                2.16*pow(config.digits, 2.0/3.0)*pow(wavenum*nodeLeng, 1.0/3.0)));
                
        Ls.push_back(floor(0.50*tau)); // TODO: Find optimal formula

        // Construct thetas
        const int nth = tau+1;
        auto [nodes, weights] = Interp::gaussLegendre(nth, EPS_NR, 0.0, PI);

        // Absorb sin(theta) into weights
        std::transform(weights.begin(), weights.end(), nodes.begin(), weights.begin(),
            [](double weight, double theta) { return weight * sin(theta); }
        );

        thetas.push_back(nodes);
        thetaWeights.push_back(weights);

        // Construct phis
        const int nph = 2*nth;
        realVec phis_lvl(nph);

        for (int iph = 0; iph < nph; ++iph)
            phis_lvl[iph] = 2.0*PI*iph/static_cast<double>(nph);

        phis.push_back(phis_lvl);

        std::cout << "   (" << lvl << "," << nth << "," << nph << ")\n";
    }
}

/* buildInteractionList()
 * Find interaction nodes
 */
void Node::buildInteractionList() {
    assert(!isRoot());
    assert(!nbors.empty());

    auto notContains = 
        [](const NodeVec& vec, const std::shared_ptr<Node> val) {
        return std::find(vec.begin(), vec.end(), val) == vec.end();
    };

    for (const auto& baseNbor : base->nbors) {
        if (baseNbor->isSrcless()) continue; // TODO: double check

        if (baseNbor->isNodeType<Leaf>() && notContains(nbors, baseNbor)) {
            leafIlist.push_back(baseNbor);
            continue;
        }

        for (const auto& branch : baseNbor->branches)
            if (notContains(nbors, branch) && !branch->isSrcless()) // TODO: double check
                iList.push_back(branch);

    }

    assert(iList.size() <= pow(6, DIM) - pow(3, DIM));

    // std::cout << iList.size() << '\n';
}

/* pushSelfToNearNonNbors()
 * Add this node to list 3 of leaf.
 * (if leaf is in list 4 of self, self is in list 3 of leaf) 
 */
void Node::pushSelfToNearNonNbors() {

    if (leafIlist.empty()) return;

    for (const auto& node : leafIlist) {
        auto leaf = dynamic_pointer_cast<Leaf>(node);

        leaf->pushToNearNonNbors(getSelf()); // call shared_from_this()

        nonNearPairs.emplace_back(getSelf(), node); // record list3-list4 pair
    }
}

void Node::buildNonNearRads() {

    for (const auto& pair : nonNearPairs) {
        auto [obsNode, srcNode] = pair;

        const size_t numObss = obsNode->srcs.size(), numSrcs = srcNode->srcs.size();

        auto nodePairRads = cmplxVec(numObss*numSrcs);

        int pairIdx = 0;
        for (size_t obsIdx = 0; obsIdx < numObss; ++obsIdx) {
            for (size_t srcIdx = 0; srcIdx < numSrcs; ++srcIdx) {
                const auto obs = obsNode->srcs[obsIdx], src = srcNode->srcs[srcIdx];

                nodePairRads[pairIdx++] = obs->getIntegratedRad(src);
            }
        }

        obsNode->nonNearRads.push_back(nodePairRads);
    }
}

/* buildMpoleToLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 */
void Node::buildMpoleToLocalCoeffs() {

    const auto [nth, nph] = getNumAngles(level); 
    localCoeffs.resize(nth*nph, vec2cd::Zero()); // TODO: Allocate elsewhere

    if (iList.empty()) return;

    for (const auto& node : iList) {

        const auto& mpoleCoeffs = node->coeffs;

        const auto& dX = center - node->center;

        const auto& transl_dX = tables.transl[level].at(dX/nodeLeng);

        for (int idx = 0; idx < nth*nph; ++idx)
            localCoeffs[idx] += transl_dX[idx] * mpoleCoeffs[idx];
    }
}

/* evalLeafIlistSols()
 * (S2L/S2T) Add contribution from list 4 to local coeffs
 */
void Node::evalLeafIlistSols() {
    for (const auto& node : leafIlist)
        evalPairSols(node, nonNearRads[nonNearPairIdx++]);
    return;

    /* No psi LUT
    const int nearIdx = std::floor((nps-1) * psi / PI);

    realVec psis_;
    for (int ips = nearIdx+1-order; ips <= nearIdx+order; ++ips)
        psis_.push_back(PI*ips/static_cast<double>(nps-1));

    for (int ips = nearIdx+1-order, k = 0; k < 2*order; ++ips, ++k) {
        const int ips_flipped = Math::flipIdxToRange(ips, nps);

        translCoeff +=
            transls[ips_flipped]
            * Interp::evalLagrangeBasis(psi,psis_,k);
    }
    */
}

/* (S2T) Evaluate sols at sources in this node due to sources in srcNode
 * and vice versa
 * srcNode : source node
 * rads    : precomputed radiation coefficients
 */
void Node::evalPairSols(const std::shared_ptr<Node> srcNode, const cmplxVec& rads) {

    const int numObss = srcs.size(), numSrcs = srcNode->srcs.size();

    cmplxVec solAtObss(numObss, 0.0);
    cmplxVec solAtSrcs(numSrcs, 0.0);

    int pairIdx = 0;
    for (size_t obsIdx = 0; obsIdx < numObss; ++obsIdx) {
        for (size_t srcIdx = 0; srcIdx < numSrcs; ++srcIdx) {
            const auto obs = srcs[obsIdx], src = srcNode->srcs[srcIdx];

            const cmplx rad = rads[pairIdx++];

            solAtObss[obsIdx] += (*currents)[src->getIdx()] * rad;
            solAtSrcs[srcIdx] += (*currents)[obs->getIdx()] * rad;
        }
    }

    for (int n = 0; n < numObss; ++n)
        (*sols)[srcs[n]->getIdx()] += Phys::C * wavenum * solAtObss[n];

    for (int n = 0; n < numSrcs; ++n)
        (*sols)[srcNode->srcs[n]->getIdx()] += Phys::C * wavenum * solAtSrcs[n];
}

// evalSelfSols without precomputed rads
void Node::evalSelfSolsDir() {

    const int numSrcs = srcs.size();

    cmplxVec solAtObss(numSrcs, 0.0);

    // TODO: Handle self-interactions
    int pairIdx = 0;
    for (size_t obsIdx = 1; obsIdx < numSrcs; ++obsIdx) { // obsIdx = 0
        for (size_t srcIdx = 0; srcIdx < obsIdx; ++srcIdx) { // srcIdx <= obsIdx 
            auto obs = srcs[obsIdx], src = srcs[srcIdx];

            const auto glObsIdx = obs->getIdx(), glSrcIdx = src->getIdx();

            const cmplx rad = obs->getIntegratedRad(src);

            solAtObss[obsIdx] += (*currents)[src->getIdx()] * rad;
            solAtObss[srcIdx] += (*currents)[obs->getIdx()] * rad;
        }
    }

    for (int n = 0; n < numSrcs; ++n)
        (*sols)[srcs[n]->getIdx()] += Phys::C * wavenum * solAtObss[n];
}