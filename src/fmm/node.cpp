#include "node.h"

/* Node(particles,branchIdx,base)
 * particles : list of particles contained in this node
 * branchidx : index of this node relative to its base node
 * base      : pointer to base node
 */
FMM::Node::Node(
    const SrcVec& srcs,
    const int branchIdx,
    Node* const base)
    : srcs(srcs), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? config.rootLeng : base->nodeLeng/2.0),
    level(base == nullptr ? 0 : base->level + 1),
    center(base == nullptr ? zeroVec :
        base->center + nodeLeng/2.0 * Math::idx2pm(branchIdx))
{
    ++numNodes;
}

/* buildInteractionList()
 * Find interaction nodes
 */
void FMM::Node::buildInteractionList() {
    assert(!isRoot());
    assert(!nbors.empty());

    auto notContains = 
        [](const NodeVec& vec, const std::shared_ptr<Node> val) {
        return std::find(vec.begin(), vec.end(), val) == vec.end();
    };

    for (const auto& baseNbor : base->nbors) {
        if (baseNbor->isSrcless()) continue;

        if (baseNbor->isNodeType<Leaf>() && notContains(nbors, baseNbor)) {
            leafIlist.push_back(baseNbor);
            continue;
        }

        for (const auto& branch : baseNbor->branches)
            if (notContains(nbors, branch) && !branch->isSrcless())
                iList.push_back(branch);
    }

    assert(iList.size() <= pow(6, DIM) - pow(3, DIM));
}

/* pushSelfToNearNonNbors()
 * Add this node to list 3 of leaf.
 * (if leaf is in list 4 of self, self is in list 3 of leaf) 
 */
void FMM::Node::pushSelfToNearNonNbors() {
    if (leafIlist.empty()) return;

    for (const auto& node : leafIlist) {
        auto leaf = dynamic_pointer_cast<Leaf>(node);

        leaf->pushToNearNonNbors(getSelf());
        nonNearPairs.emplace_back(leaf, getSelf()); // record list4-list3 pair
    }
}

/* translateCoeffs()
 * (M2L) Translate mpole coeffs of interaction nodes into local coeffs at center,
 * then apply integration weights for anterpolation
 */
void FMM::Node::translateCoeffs() {
    if (iList.empty()) return;

    localCoeffs.fillZero();

    const auto& transl = tables[level].transl;
    const size_t nDir = localCoeffs.size();

    for (const auto& node : iList) {
        const auto& dX = center - node->center;
        const auto& transl_dX = transl.at(dX/nodeLeng);

        Eigen::Map<arrXcd> mpoleTheta(node->coeffs.theta.data(), nDir);
        Eigen::Map<arrXcd> mpolePhi(node->coeffs.phi.data(), nDir);

        Eigen::Map<arrXcd> thetaLocal(localCoeffs.theta.data(), nDir);
        Eigen::Map<arrXcd> phiLocal(localCoeffs.phi.data(), nDir);
    
        thetaLocal += transl_dX * mpoleTheta;
        phiLocal += transl_dX * mpolePhi;
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

/* evalLeafIlistSols()
 * (S2L/S2T) Add contribution from list 4 to local coeffs
 */
void FMM::Node::evalLeafIlistSols() {
    // Do nothing! Contribution from list 4 node is 
    // to be evaluated by Leaf::evalNearNonNborSols()
    //for (const auto& node : leafIlist)
    //    evalPairSols(node, nonNearRads[nonNearPairIdx++]);
    //return;
}

void FMM::Node::printFarFld(const std::string& fname) {
    namespace fs = std::filesystem;
    fs::path dir = "out/ff";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream farfile(dir/fname);
    farfile << std::setprecision(15) << std::scientific;

    const auto& angles_lvl = angles[level];
    size_t nDir = angles_lvl.getNumAllAngles();

    for (int iDir = 0; iDir < nDir; ++iDir) {
        const auto& krhat = angles_lvl.khat[iDir] * k;

        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += (*currents)[src->getIdx()] * src->getFarAlongDir(krhat);

        const vec3cd& far = Phys::C * k * angles_lvl.ImRR[iDir] * dirFar;

        farfile << far << '\n';
    }

    // Also print out angles (coordinates of farsols)
    std::ofstream thfile(dir/"thetas.txt"), phfile(dir/"phis.txt");
    angles[level].printAngles(thfile, phfile);
}