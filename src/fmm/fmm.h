#pragma once

#include "clock.h"
#include "types.h"

extern ClockTimes t;

namespace FMM {
    constexpr int DIM = 3;
    constexpr int numDir = 26;

    enum class Dir {
        W, E, S, N, D, U,
        SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
        DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
    };

    // TODO: Move to angles.h
    struct Angles {
        Angles() = default;

        pair2i getNumAngles(const int level) const {
            return std::make_pair(thetas[level].size(), phis[level].size());
        }

        std::vector<realVec> thetas;
        std::vector<realVec> thetaWeights;
        std::vector<realVec> phis;
        std::vector<int> Ls;
    };


    class Node;
    class Stem;
    class Leaf;

    class Tables;
}