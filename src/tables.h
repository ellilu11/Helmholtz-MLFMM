#pragma once

#include "interp.h"
#include "map.h"

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
    std::vector<Map<vecXcd>> transl;
    std::vector<HashMap<vecXcd>> interpPsi;
    std::vector<HashMap<int>> ssps;

};

