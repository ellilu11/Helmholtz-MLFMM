#pragma once

#include "../interp.h"
#include "../maps.h"
#include "fmm.h"

class FMM::Tables {

public:
    Tables() = default;

    Tables(const Angles& angles, const Config& config, double wavenum, int maxLevel)
        : angles(angles), 
          wavenum(wavenum), 
          maxLevel(maxLevel),
          order(config.interpOrder), 
          rootLeng(config.rootLeng), 
          overInterp(config.overInterp),
          dists(Math::getINodeDistances()),
          rhats(Math::getINodeDirections())
    {
        buildAngularTables();

        buildInterpTables();

        buildTranslationTable();
    }

    // Angular tables
    std::vector<std::vector<vec3d>> khat;
    std::vector<std::vector<mat23d>> toThPh;
    std::vector<std::vector<mat3d>> ImRR;

    // M2M interpolation tables
    std::vector<std::vector<interpPair>> interpTheta;
    std::vector<std::vector<interpPair>> interpPhi;

    // L2L interpolation tables
    std::vector<std::vector<interpPair>> invInterpTheta;
    std::vector<std::vector<interpPair>> invInterpPhi;

    // M2L translation table
    std::vector<VecHashMap<vecXcd>> transl;

private:
    std::vector<interpPair> getInterpThetaAtLvl(int, int);

    std::vector<interpPair> getInterpPhiAtLvl(int, int);

    Map<vecXcd> getAlphaAtLvl(int);

    HashMap<interpPair> getInterpPsiAtLvl(int);

    void buildAngularTables();

    void buildInterpTables();

    void buildTranslationTable();

    Angles angles;
    double wavenum;
    int maxLevel;
    int order;
    double rootLeng;
    double overInterp;

    realVec dists;
    std::vector<vec3d> rhats;
};
