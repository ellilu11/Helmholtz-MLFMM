#pragma once

#include "node.h"

struct FMM::NearPair {
    friend class Nearfield;

    NearPair(
        std::shared_ptr<Node> obsLeaf, 
        std::shared_ptr<Node> srcLeaf)
        : pair(std::move(obsLeaf), std::move(srcLeaf)) 
    {}

    matXcd getNearMatrix() const;

    NodePair pair;
    std::vector<cmplx> rads;
};

class FMM::Nearfield {

public:
    Nearfield(int);

    void evaluateSols();

    // sparseMat<cmplx> getNearMatrix(int) const;

    std::vector<NearPair> getNearPairs() const { return nearPairs; }

    std::vector<NearPair> getSelfPairs() const { return selfPairs; }

private: 
    void findNodePairs();

    void buildTriPairs();

    void buildPairRads(std::vector<Eigen::Triplet<cmplx>>&);

    void buildSelfRads(std::vector<Eigen::Triplet<cmplx>>&);

    static void evalPairSols(const NearPair&);

    static void evalSelfSols(const NearPair&);

public:
    sparseMat<cmplx> nearMat;

private:
    std::vector<NearPair> nearPairs;
    std::vector<NearPair> selfPairs;
};