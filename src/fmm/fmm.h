#pragma once

#include "../clock.h"
#include "../config.h"
#include "../excite.h"
#include "../interp.h"
#include "../phys.h"
#include "../types.h"

extern ClockTimes t;

namespace FMM {
    // Constants
    constexpr int DIM = 3;
    constexpr int numDir = 26;
    constexpr double EPS_NR = 1.0E-9; // Newton-Raphson precision

    // Types
    enum class Dir {
        W, E, S, N, D, U,
        SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
        DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
    };

    class Node;
    class Stem;
    class Leaf;
    struct Angles;
    class Tables;

    using NodeVec = std::vector<std::shared_ptr<Node>>;
    using NodePair = std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>>;

    using LeafVec = std::vector<std::shared_ptr<Leaf>>;
    using LeafPair = std::pair<std::shared_ptr<Leaf>, std::shared_ptr<Leaf>>;

    // Global data
    Config config;
    double wavenum;
    std::vector<Angles> angles;
    std::vector<Tables> tables;
    inline int maxLevel = 0;

    LeafVec leaves;
    std::vector<NodePair> nonNearPairs;
    std::vector<LeafPair> nearPairs;

    // Move to global namespace?
    std::shared_ptr<vecXcd> lvec;
    std::shared_ptr<vecXcd> rvec;
    std::shared_ptr<vecXcd> currents;

    // Functions
    void initGlobal(
        const Config&,
        const std::shared_ptr<Excitation::PlaneWave>&,
        int);

    void buildTables();

    void resetLeaves() {
        leaves.clear();
        nearPairs.clear();
    }
}