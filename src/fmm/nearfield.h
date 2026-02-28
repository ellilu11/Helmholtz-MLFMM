#pragma once

#include "node.h"

struct FMM::NearPair {
    friend class Nearfield;

    NearPair(
        std::shared_ptr<Node> obsLeaf, 
        std::shared_ptr<Node> srcLeaf)
        : pair(std::move(obsLeaf), std::move(srcLeaf)) 
    {}

    NodePair pair;
    std::vector<cmplx> efie;
    std::vector<pair2cd> mfie;
};

class FMM::Nearfield {

public :
    Nearfield();

    void evaluateSols();

private : 
    void findPairs();

    void buildPairRads();

    void buildSelfRads();

    static void evalPairSols(const NearPair&);

    static void evalSelfSols(const NearPair&);

private:
    std::vector<NearPair> nearPairs;
    std::vector<NearPair> selfPairs;
};