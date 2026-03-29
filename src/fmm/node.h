#pragma once

#include <cassert>
#include <numeric>
#include <queue>

#include "coeffs.h"
#include "fmm.h"
#include "level.h"
#include "../source.h"

class FMM::Node : public std::enable_shared_from_this<Node> {
    friend class Farfield;
    friend class Nearfield;

public:
    Node(const SrcVec&, const int, Node* const, bool);

    void postProcess();

    bool isRoot() const { return base == nullptr; }

    bool isLeaf() const { return branches.empty(); }
    
    bool isSrcless() const { return srcs.empty(); }

private:    
    void subdivide();

    void buildNeighbors();

    void buildInteractionList();
    
    void pushSelfToNearNonNbors();

    void pushToNearNonNbors(const std::shared_ptr<Node>& node) {
        nearNonNbors.push_back(node);
    }

    std::shared_ptr<Node> getNeighborGeqSize(const Dir) const;

    NodeVec getNeighborsLeqSize(const std::shared_ptr<Node>, const Dir) const;

    void findTris();

private:
    std::vector<Coeffs> radPats;
    std::vector<Coeffs> recPatsH;

    Coeffs coeffs;
    Coeffs localCoeffs;

    NodeVec branches;
    NodeVec nbors;
    NodeVec iList; // list 2
    NodeVec nearNbors; // list 1
    NodeVec nearNonNbors; // list 3
    NodeVec leafIlist; // list 4

    std::vector<int> iTris; // indices of unique triangles of RWGs in this node

    SrcVec srcs;
    const int branchIdx;
    Node* const base;
    const double nodeLeng; // declare before center!
    const vec3d center;
    const int lvl;
};