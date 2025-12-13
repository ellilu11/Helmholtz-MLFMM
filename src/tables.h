#pragma once

#include <map>
#include "interp.h"

constexpr double DistEPS = 1.0E-3;
constexpr double PsiEPS = 1.0E-9;

struct DistComp {
    bool operator()(double x, double y) const noexcept {
        return x + DistEPS < y;
    }
};

struct PsiHash {
    std::size_t operator()(double x) const noexcept {
        long long q = static_cast<long long>(std::llround(x / PsiEPS));
        return std::hash<long long>{}(q);
    }
};

struct PsiEq {
    bool operator()(double x, double y) const noexcept {
        return std::fabs(x - y) <= PsiEPS;
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

    std::vector<std::vector<mat23d>> matToThPh;
    std::vector<std::vector<mat32d>> matFromThPh;
    
    // std::vector<std::vector<mat3d>> matToSph;
    // std::vector<std::vector<mat3d>> matFromSph;

    // M2M and L2L interpolation tables
    std::vector<std::vector<realVec>> interpTheta;
    std::vector<std::vector<int>> ts;

    std::vector<std::vector<realVec>> interpPhi;
    std::vector<std::vector<int>> ss;

    // M2L translation tables
    std::vector<std::map<double,vecXcd,DistComp>> transl;
    std::vector<std::unordered_map<double,vecXd,PsiHash,PsiEq>> interpPsi;
    std::vector<std::unordered_map<double,int,PsiHash,PsiEq>> ssps;

};

