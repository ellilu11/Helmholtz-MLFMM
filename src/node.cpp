#include "node.h"

Config Node::config;
double Node::wavenum;
std::vector<realVec> Node::thetas;
std::vector<realVec> Node::thetaWeights;
std::vector<realVec> Node::phis;
// std::vector<realVec> Node::psis;
Tables Node::tables;

NodeVec Node::nodes;

void Node::setNodeParams(
    const Config& config_, const std::shared_ptr<Source>& Einc) {

    config = config_;
    wavenum = Einc->wavenum;
}

/* Node(particles,branchIdx,base)
 * particles : list of particles contained in this node
 * branchidx : index of this node relative to its base node
 * base      : pointer to base node
 */
Node::Node(
    const RWGVec& rwgs,
    const int branchIdx,
    Node* const base)
    : rwgs(rwgs), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? config.rootLeng : base->nodeLeng / 2.0),
    level(base == nullptr ? 0 : base->level + 1),
    center(base == nullptr ? zeroVec :
        base->center + nodeLeng / 2.0 * Math::idx2pm(branchIdx)),
    label(0)
{
    nodeIdx = numNodes++;
}

/* buildAngularSamples()
 * Compute theta and phi samples at each level
 */
void Node::buildAngularSamples() {

    constexpr double EPS = 1.0E-9;

    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        const double nodeLeng = config.rootLeng / pow(2.0, lvl);

        // Use excess bandwidth formula
        const int L = ceil(
                (1.73*wavenum*nodeLeng +
                2.16*pow(config.digits, 2.0/3.0)*pow(wavenum*nodeLeng, 1.0/3.0)));

        // Construct thetas
        const int nth = L+1;
        const auto [nodes, weights] = Interp::gaussLegendre(nth, EPS, 0.0, PI);

        thetas.push_back(nodes);
        thetaWeights.push_back(weights);

        // Construct phis
        const int nph = 2*(L+1);
        realVec phis_lvl;

        for (int iph = 0; iph < nph; ++iph)
            phis_lvl.push_back(2.0*PI*iph/static_cast<double>(nph));

        phis.push_back(phis_lvl);

        // Construct psis
        /*const int nps = std::floor(Q*L);
        realVec psis_lvl;

        for (int ips = 0; ips < nps; ++ips)
            psis_lvl.push_back(PI*ips/static_cast<double>(nps-1));

        psis.push_back(psis_lvl);*/

        std::cout << "   (Lvl,Nth,Nph) = "
                  << "(" << lvl << "," << nth << "," << nph << ")\n";
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

    // std::cout << iList.size() << ' ';
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
    }
}

/* buildMpoleToLocalCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center
 */
void Node::buildMpoleToLocalCoeffs() {
    if (isRoot()) return;

    const int order = config.interpOrder;

    const auto [nth, nph] = getNumAngles(level);
    localCoeffs.resize(nth*nph, vec3cd::Zero());

    const auto nps = std::floor(Q*(nth-1));

    // auto start = Clock::now();

    for (const auto& node : iList) {
        const auto& mpoleCoeffs = node->coeffs;

        const auto& dR = center - node->center;
        const auto r = dR.norm();
        const auto& rhat = dR / r;

        assert(nodeLeng == node->nodeLeng);
        const auto& translVec = tables.transl[level].at(r / nodeLeng);

        size_t idx = 0;
        for (int ith = 0; ith < nth; ++ith) {

            for (int iph = 0; iph < nph; ++iph) {

                const auto& khat = tables.khat[level][idx];

                const double psi = acos(khat.dot(rhat));
                // const double psi = khat.dot(rhat);

                // start = Clock::now();

                const int s = tables.ssps[level].at(psi);

                cmplx translCoeff = 0.0;

                /* No psi LUT
                realVec psis_;
                for (int ips = s+1-order; ips <= s+order; ++ips)
                    psis_.push_back(PI*ips/static_cast<double>(nps-1));

                for (int ips = s+1-order, k = 0; k < 2*order; ++ips, ++k) {
                    const int ips_flipped = Math::flipIdxToRange(ips, nps);

                    translCoeff +=
                        translVec[ips_flipped]
                        * Interp::evalLagrangeBasis(psi,psis_,k);
                }*/
                
                const auto& interpVec = tables.interpPsi[level].at(psi);

                // t.M2L_lookup += Clock::now() - start;

                for (int ips = s+1-order, k = 0; k < 2*order; ++ips, ++k) {

                    // using cos(-psi) = cos(psi), cos(2pi-psi) = cos(psi)
                    const int ips_flipped = Math::flipIdxToRange(ips, nps); 

                    translCoeff +=
                         translVec[ips_flipped] * interpVec[k];
                }

                localCoeffs[idx] += translCoeff * mpoleCoeffs[idx];

                idx++;
            }
        }
    }
}

/* evalLeafIlistSols()
 * (S2L) Add contribution from list 4 to local coeffs
 */
void Node::evalLeafIlistSols() {

}

/* evalPairSols(srcNode)
 * (S2T) Evaluate sols at RWGs in this node due to RWGs in srcNode
 * srcNode : source node
 */
void Node::evalPairSols(const std::shared_ptr<Node> srcNode) {

    const auto& srcs = srcNode->rwgs;

    for (const auto& obs : rwgs) {
        // Integrate E_rad over obs RWG
        cmplx sol = 0;

        for (const auto& src : srcs)
            // Integrate G \dot J over src RWG
             sol += obs->getIntegratedRad(src, wavenum);

        obs->addToSol(C * wavenum * sol);
    }
}

/* evalSelfSols()
 * (S2T) Evaluate sols at RWGs in this node due to other RWGs
 * in this node
 */
void Node::evalSelfSols() {

    for (const auto& obs : rwgs) {
        // Integrate E_rad over obs RWG
        cmplx sol = 0;

        for (const auto& src : rwgs) {
            if (src == obs) continue; // TODO: Use analytic expression

            // Integrate G \dot J over src RWG
            sol += obs->getIntegratedRad(src, wavenum);

        }

        obs->addToSol(C * wavenum * sol);
    }
}

/* getFarSol()
 * Return sols at all sampled directions at distance r 
 * due to all RWGs in this node using farfield approximation
 */
std::vector<vec3cd> Node::getFarSols(double r) {

    assert(r >= 5.0 * config.rootLeng); // verify farfield condition

    const cmplx C = -iu * c0 * wavenum * mu0
        * exp(iu*wavenum*r) / (4.0*PI*r);

    const auto [nth, nph] = getNumAngles(level);

    std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const auto& ImKK = tables.ImKK[level][idx];
            const auto& kvec = tables.khat[level][idx] * wavenum;

            vec3cd dirCoeff = vec3cd::Zero();
            for (const auto& rwg : rwgs)
                dirCoeff += rwg->getRadAlongDir(center, kvec);

            sols[idx] = C * ImKK * dirCoeff;

            idx++;
        }
    }

    return sols;
}

/* getFarSolsFromCoeffs(r)
 * Return sols at all sampled directions at distance r
 * from mpole coefficients of this node
 */
std::vector<vec3cd> Node::getFarSolsFromCoeffs(double r) {
    assert(!rwgs.empty());

    const cmplx C = -iu * c0 * wavenum * mu0
        * exp(iu*wavenum*r) / (4.0*PI*r);

    const auto [nth, nph] = getNumAngles(level);

    std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const auto kvec = tables.khat[level][idx] * wavenum;

            sols[idx] =
                C * exp(-iu*kvec.dot(center))
                // * tables.matFromSph[level][idx]
                * coeffs[idx];

            idx++;
        }
    }

    return sols;
}