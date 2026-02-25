#pragma once

#include <cassert>
#include <numeric>
#include <queue>

#include "angles.h"
#include "coeffs.h"
#include "fmm.h"
#include "tables.h"
#include "../source.h"

class FMM::Node : public std::enable_shared_from_this<Node> {
    friend class Nearfield;

public:
    Node(const SrcVec&, const int, Node* const, bool);

    void buildLists();

    void resizeCoeffs();

    static void buildRadPats();

    Coeffs buildMpoleCoeffs();

    Coeffs mergeMpoleCoeffs();

    void buildLocalCoeffs();

    static void evaluateSols();

    static void addInterpCoeffs(const Coeffs&, Coeffs&, int, int);

    static void addAnterpCoeffs(const Coeffs&, Coeffs&, int, int);

    void printScatteredField(const std::string&, int, bool = 0);

    static int getMaxLvl() { return maxLevel; }

    static int getNumNodes() { return numNodes; }

    SrcVec getSrcs() const { return srcs; }
    
    int getBranchIdx() const { return branchIdx; }

    Node* getBase() const { return base; }
    
    double getNodeLeng() const { return nodeLeng; }
    
    int getLevel() const { return level; }

    vec3d getCenter() const { return center; }
    
    Coeffs getMpoleCoeffs() const { return coeffs; }
    
    Coeffs getLocalCoeffs() const { return localCoeffs; }

    bool isRoot() const { return base == nullptr; }

    bool isLeaf() const { return branches.empty(); }
    
    bool isSrcless() const { return srcs.empty(); }

private:
    std::shared_ptr<Node> getNeighborGeqSize(const Dir) const;

    NodeVec getNeighborsLeqSize(const std::shared_ptr<Node>, const Dir) const;
    
    void subdivideNode();

    void buildNeighbors();

    void balanceNeighbors();

    void buildInteractionList();
    
    void pushSelfToNearNonNbors();

    Coeffs getShiftedLocalCoeffs(int) const;

    void translateCoeffs();

    void evalFarSols();

    void pushToNearNonNbors(const std::shared_ptr<Node>& node) {
        nearNonNbors.push_back(node);
    }

private:
    inline static int numNodes = 0;

    std::vector<Coeffs> radPats;
    Coeffs coeffs;
    Coeffs localCoeffs;

    NodeVec branches;
    NodeVec nbors;
    NodeVec iList; // list 2
    NodeVec nearNbors; // list 1
    NodeVec nearNonNbors; // list 3
    NodeVec leafIlist; // list 4

    SrcVec srcs;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const int level;
    const vec3d center;
};