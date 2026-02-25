#include "node.h"

/* Node(particles,branchIdx,base)
 * particles : list of particles contained in this node
 * branchidx : index of this node relative to its base node
 * base      : pointer to base node
 */
FMM::Node::Node(
    const SrcVec& srcs,
    const int branchIdx,
    Node* const base,
    bool buildLeaf)
    : srcs(srcs), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? config.rootLeng : base->nodeLeng/2.0),
    level(base == nullptr ? 0 : base->level + 1),
    center(base == nullptr ? zeroVec :
        base->center + nodeLeng/2.0 * Math::idx2pm(branchIdx))
{
    if (buildLeaf) maxLevel = std::max(level, maxLevel);
    else subdivideNode();

    ++numNodes;
}

void FMM::Node::subdivideNode() {
    // Assign every src to a branch based on src center relative to node center
    std::array<SrcVec, 8> branchSrcs;
    for (const auto& src : srcs) {
        size_t idx = Math::bools2Idx(src->getCenter() > center);
        branchSrcs[idx].push_back(src);
    }

    // Construct branch nodes
    branches.resize(8);
    for (size_t k = 0; k < 8; ++k) {
        bool buildLeaf = branchSrcs[k].size() <= config.maxNodeSrcs;
        branches[k] = std::make_shared<Node>(branchSrcs[k], k, this, buildLeaf);
    }
}

/* buildNeighbors()
 * Find all neighbor nodes of equal or greater size
 * If node is leaf, find all neighbor leaves of equal or lesser size (list 1)
 */
void FMM::Node::buildNeighbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);
        if (!nbor) continue; // continue if nbor is nullptr
        
        nbors.push_back(nbor);

        if (!isLeaf()) continue;
        auto nbors = getNeighborsLeqSize(nbor, dir);
        nearNbors.insert(nearNbors.end(), nbors.begin(), nbors.end());
    }

    assert(nbors.size() <= numDir);
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

        if (baseNbor->isLeaf() && notContains(nbors, baseNbor)) {
            leafIlist.push_back(baseNbor); // TODO: Subdivide baseNbor instead of adding to leafIlist
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
        node->pushToNearNonNbors(shared_from_this());
        nonNearPairs.emplace_back(node, shared_from_this()); // record list4-list3 pair
    }
}

/* buildLists()
 * Find neighbor and interaction lists.
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void FMM::Node::buildLists() {
    if (!isRoot()) {
        buildNeighbors();
        buildInteractionList();
        pushSelfToNearNonNbors();
    }

    if (isLeaf()) leaves.push_back(shared_from_this());

    for (const auto& branch : branches)
        branch->buildLists();
}

void FMM::Node::resizeCoeffs() {
    const auto [nth, nph] = angles[level].getNumAngles();

    coeffs.resize(nth*nph);
    localCoeffs.resize(nth*nph);

    for (const auto& branch : branches)
        branch->resizeCoeffs();
}