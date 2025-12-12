#pragma once

#include <map>
#include "interp.h"

constexpr double DistEPS = 1.0E-3;
constexpr double PsiEPS = 1.0E-9;

struct DistComp {
    bool operator()(double x, double y) const {
        return x + DistEPS < y;
    }
};

struct PsiComp {
    bool operator()(double x, double y) const {
        return x + PsiEPS < y;
    }
};

struct Tables {
    Tables() = default;
    Tables(const Config& config)
    {
        order = config.interpOrder;

        buildAngularTables();
        
        buildInterpThetaTable();
        buildInterpPhiTable();
        
        buildTranslationTable();
        buildInterpPsiTable();
    }
    
    void buildAngularTables();

    void buildInterpThetaTable();
    void buildInterpPhiTable();

    void buildTranslationTable();
    void buildInterpPsiTable();

    int order;

    // Angular tables
    std::vector<std::vector<mat3d>> ImKK;
    std::vector<std::vector<vec3d>> khat;

    std::vector<std::vector<Eigen::Matrix<double,2,3>>> matToThPh;
    std::vector<std::vector<Eigen::Matrix<double,3,2>>> matFromThPh;
    
    std::vector<std::vector<mat3d>> matToSph;
    std::vector<std::vector<mat3d>> matFromSph;

    // M2M and L2L interpolation tables
    std::vector<std::vector<realVec>> interpTheta;
    std::vector<std::vector<int>> ts;

    std::vector<std::vector<realVec>> interpPhi;
    std::vector<std::vector<int>> ss;

    // M2L translation tables
    std::vector<std::map<double,cmplxVec,DistComp>> transl;
    std::vector<std::map<double,vecXd,PsiComp>> interpPsi;
    // std::vector<std::map<double,int>> ssps;

};

