#pragma once

#include "math.h"

struct Tables {
    Tables() = default;
    Tables(const int maxLevel,
        const int order,
        const double wavenum, 
        std::vector<realVec>& thetas, 
        std::vector<realVec>& phis) 
    {
        buildAngularTables(maxLevel, wavenum, thetas, phis);
        buildInterpThetaTable(maxLevel, order, thetas);
        buildInterpPhiTable(maxLevel, order, phis);
    }

    void buildAngularTables(
        const int, const double, const std::vector<realVec>&, const std::vector<realVec>&);

    void buildInterpThetaTable(const int, const int, const std::vector<realVec>&);

    void buildInterpPhiTable(const int, const int, const std::vector<realVec>&);

    // void buildYlmTables(const int);

    // Angular tables
    std::vector<std::vector<mat3d>> ImKK;
    std::vector<std::vector<vec3d>> kvec;
    std::vector<std::vector<mat23d>> matToThPh;

    // Lagrange interpolation tables
    std::vector<std::vector<realVec>> interpTheta;
    std::vector<std::vector<size_t>> T;

    std::vector<std::vector<realVec>> interpPhi;
    std::vector<std::vector<size_t>> S;

    // Ylm tables
    /*std::vector<realVec> coeffYlm_;
    std::vector<realVec> fallingFact_;
    std::vector<realVec> legendreSum_;
    std::vector<realVec> fracCoeffYlm_;
    std::vector<realVec> A_;
    std::vector<realVec> Aexp_;*/
};
