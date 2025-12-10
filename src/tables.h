#pragma once

#include "interp.h"

struct Tables {
    Tables() = default;
    Tables(const int maxLevel,
        const int order,
        std::vector<realVec>& thetas, 
        std::vector<realVec>& phis,
        std::vector<int>& Ls,
        const double wavenum,
        const double rootLeng)
    {
        buildAngularTables(maxLevel, thetas, phis, wavenum);
        
        buildInterpThetaTable(maxLevel, order, thetas);
        buildInterpPhiTable(maxLevel, order, phis);
        
        // buildTranslationTable(maxLevel, order, Ls, wavenum, rootLeng);
        // buildInterpPsiTable(maxLevel, order, thetas, phis);

    }

    void buildAngularTables(
        const int, const std::vector<realVec>&, const std::vector<realVec>&, const double);

    void buildInterpThetaTable(const int, const int, const std::vector<realVec>&);

    void buildInterpPhiTable(const int, const int, const std::vector<realVec>&);

    void buildTranslationTable(
        const int, const int, const std::vector<int>&, const double, const double);

    void buildInterpPsiTable(
        int, int, const std::vector<realVec>&, const std::vector<realVec>&);

    // Angular tables
    std::vector<std::vector<mat3d>> ImKK;
    std::vector<std::vector<vec3d>> kvec;
    std::vector<std::vector<Eigen::Matrix<double,2,3>>> matToThPh;
    std::vector<std::vector<Eigen::Matrix<double,3,2>>> matFromThPh;
    
    std::vector<std::vector<mat3d>> matToSph;
    std::vector<std::vector<mat3d>> matFromSph;

    // Lagrange interpolation tables
    std::vector<std::vector<realVec>> interpTheta;
    std::vector<std::vector<size_t>> ts;

    std::vector<std::vector<realVec>> interpPhi;
    std::vector<std::vector<size_t>> ss;

    std::vector<realVec> interpPsi;

    // M2L translation tables
    std::vector<std::vector<cmplxVec>> transl;
    realVec iNodeDists;
    std::vector<vec3d> iNodeDirs;

};
