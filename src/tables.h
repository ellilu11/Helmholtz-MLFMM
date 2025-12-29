#pragma once

#include "clock.h"
#include "interp.h"
#include "map.h"

struct Tables {
    Tables() = default;
    Tables(const Config& config)
        : order(config.interpOrder)
    {
        dists = Math::getINodeDistances();
        rhats = Math::getINodeDirections();

        buildAngularTables();
        
        buildInterpTables();
        
        buildTranslationTable();
    }
    
    void buildAngularTables();

    std::vector<interpPair> getInterpThetaAtLvl(int, int);

    std::vector<interpPair> getInterpPhiAtLvl(int, int);

    void buildInterpTables();

    Map<vecXcd> getAlphaAtLvl(int);

    HashMap<interpPair> getInterpPsiAtLvl(int);

    void buildTranslationTable();

    int order;
    realVec dists;
    std::array<vec3d,316> rhats;

    // Angular tables
    std::vector<std::vector<vec3d>> khat;
    std::vector<std::vector<mat23d>> toThPh;

    // M2M interpolation tables
    std::vector<std::vector<interpPair>> interpTheta;
    std::vector<std::vector<interpPair>> interpPhi;

    // L2L interpolation tables
    std::vector<std::vector<interpPair>> invInterpTheta;
    std::vector<std::vector<interpPair>> invInterpPhi;

    // M2L translation table
    std::vector<VecHashMap<vecXcd>> transl;

};

