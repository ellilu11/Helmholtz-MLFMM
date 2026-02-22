#pragma once

#include "fmm.h"

struct FMM::Angles {
    Angles() = default;

    Angles(int level) {
        buildAngularSamples(level);
        buildAngularMatrices();
    }

    void buildAngularSamples(int);

    void buildAngularMatrices();

    void printAngles(std::ofstream&, std::ofstream&);

    pair2i getNumAngles() const {
        return std::make_pair(thetas.size(), phis.size());
    }

    size_t getNumDirs() const {
        return thetas.size() * phis.size();
    }

    std::vector<vec3d> khat;
    std::vector<mat23d> toThPh;
    std::vector<mat3d> ImRR;

    std::vector<double> thetas;  // theta samples
    std::vector<double> weights; // weights of theta samples
    std::vector<double> phis;    // phi samples
    int L;           // M2L series truncation number
};