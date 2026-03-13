#include "node.h"
#include "../mesh/rwg.h"

std::vector<int> FMM::Node::numNodesPerLvl = std::vector<int>(10, 0);

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
    nodeLeng(base == nullptr ? Mesh::rootLeng : base->nodeLeng/2.0),
    level(base == nullptr ? 0 : base->level + 1),
    center(base == nullptr ? Mesh::rootCenter :
        base->center + nodeLeng/2.0 * Math::idx2pm(branchIdx))
{
    if (isRoot()) std::cout << " Building FMM tree...\n";

    if (buildLeaf) {
        // Assign contiguous global indices to srcs in this leaf node
        // TODO: Use Morton ordering to improve spatial locality of srcs in leaf nodes
        // for (const auto& src : srcs) src->setIdx(glSrcIdx++);

        // Update maxLevel if needed
        maxLevel = std::max(level, maxLevel);
    }
    else subdivideNode();

    ++numNodes;

    // For debugging: count number of nodes at each level
    assert(level < numNodesPerLvl.size());
    ++numNodesPerLvl[level];
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

/* buildNeighborsGeqSize()
 * Find all neighbor nodes of equal or greater size
 */
void FMM::Node::buildNeighbors(bool isBalanced) {
    nbors.clear();

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);
        if (!nbor) continue; // continue if nbor is nullptr

        if (isBalanced) assert(std::abs(nbor->level-level) <= 1); //
        nbors.push_back(nbor);

        if (!isLeaf() || !isBalanced) continue;
        auto nbors = getNeighborsLeqSize(nbor, dir);
        nearNbors.insert(nearNbors.end(), nbors.begin(), nbors.end());
    }
    assert(nbors.size() <= numDir);
}

void FMM::Node::doPreBalance() {
    if (isLeaf()) leaves.push_back(shared_from_this());

    if (!isRoot()) buildNeighbors(false);

    for (const auto& branch : branches)
        branch->doPreBalance();
}

void FMM::Node::balanceNodes() {
    std::queue<std::shared_ptr<Node>> queue;
    assert(!leaves.empty()); // check leaves are populated by doPreBalance()
    for (const auto& leaf : leaves) queue.push(leaf);

    while (!queue.empty()) {
        auto node = queue.front();
        assert(node->isLeaf()); // check only leaf nodes are added to queue
        queue.pop();
        
        node->buildNeighbors(false); // Rebuild neighbors

        for (const auto& nbor : node->nbors) {
            // Only subdivide nbor if it is a leaf and more than one level coarser than this node
            if (!nbor->isLeaf() || nbor->level >= node->level-1) continue;

            // Subdivide nbor and add its branches to queue
            nbor->subdivideNode();
            for (const auto& branch : nbor->branches)
                queue.push(branch);

            // Recheck this node after subdividing nbor
            queue.push(node); 
            break;
        }
    }

    leaves.clear(); // reset leaves so they can be repopulated by doPostBalance()
}

/* buildInteractionList()
 * Find interaction nodes
 */
void FMM::Node::buildInteractionList() {
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

/* doPostBalance()
 * Find neighbor and interaction lists.
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void FMM::Node::doPostBalance() {
    if (!isRoot()) {
        buildNeighbors(true);
        buildInteractionList();
        pushSelfToNearNonNbors();
    }

    if (isLeaf()) leaves.push_back(shared_from_this());

    findTris(); // TODO: Only call for leaves and non-near stems

    for (const auto& branch : branches)
        branch->doPostBalance();

    if (isRoot()) {
        std::cout << "   # Nodes: " << numNodes << '\n';
        std::cout << "   # Leaves: " << leaves.size() << '\n';
        std::cout << "   Max node level: " << maxLevel << "\n\n";

        //for (int lvl = 0; lvl <= maxLevel; ++lvl)
        //    std::cout << "   # Nodes at level " << lvl << ": " << numNodesPerLvl[lvl] << '\n';
    }
}

void FMM::Node::resizeCoeffs() {
    if (isRoot() && isLeaf()) return;

    const auto [nth, nph] = angles[level].getNumAngles();
    coeffs.resize(nth*nph);
    localCoeffs.resize(nth*nph);

    for (const auto& branch : branches)
        branch->resizeCoeffs();
}

/* findTris()
 * Find all unique triangles of RWGs in this node
 */
void FMM::Node::findTris() {
    if (isSrcless()) return;

    for (const auto& src : srcs) {
        if (!src->isSrcType<Mesh::RWG>()) continue;
        const auto rwg = dynamic_pointer_cast<Mesh::RWG>(src);

        auto [iTri0, iTri1] = rwg->getTrisIdx();
        iTris.push_back(iTri0);
        iTris.push_back(iTri1);
    }

    std::sort(iTris.begin(), iTris.end());
    iTris.erase(std::unique(iTris.begin(), iTris.end()), iTris.end());
}
