#pragma once

#include "../maps.h"
#include "../math.h"
#include "fmm.h"

struct FMM::Tables {

public:
    Tables() = default;

    Tables(int level) : level(level)
    {
        if (level < maxLevel) buildInterpTables();
        buildTranslationTable();
    }

    static void buildDists() {
        dists = Math::getINodeDistances();
        rhats = Math::getINodeDirections();
        dXs = Math::getINodeDistVecs();
    }

    static void clearDists() {
        dists.clear();
        rhats.clear();
        dXs.clear();
    }

private:
    std::vector<interpPair> getInterpTheta(int, int);

    std::vector<interpPair> getInterpPhi(int, int);

    Map<std::vector<cmplx>> getAlpha();

    HashMap<interpPair> getInterpPsi();

    void buildInterpTables();

    void buildTranslationTable();

public:
    // M2M interpolation tables
    std::vector<interpPair> interpTheta;
    std::vector<interpPair> interpPhi;

    // M2L translation table
    VecHashMap<arrXcd> transl;

private:
    static std::vector<double> dists;
    static std::vector<vec3d> rhats;
    static std::vector<vec3d> dXs;

    int level;
};
