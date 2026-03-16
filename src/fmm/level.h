#pragma once

#include "../maps.h"
#include "../math.h"
#include "fmm.h"

class FMM::Level {

public:
    Level() = default;

    Level(int level)
        : level(level)
    {
        buildAngularSamples();
        buildAngularMatrices();
    }

    void buildAngularSamples();

    void buildAngularMatrices();

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

    void buildTables(int maxLevel) {
        if (level < maxLevel) {
            interpTheta = getInterpTheta(level+1, level);
            interpPhi = getInterpPhi(level+1, level);
        }
            
        buildTranslationTable();
    }

    pair2i getNumAngles() const {
        return std::make_pair(thetas.size(), phis.size());
    }

    size_t getNumDirs() const {
        return thetas.size() * phis.size();
    }

    void printAngles(std::ofstream& thfile, std::ofstream& phfile) {
        thfile << std::setprecision(15);
        for (const auto& theta : thetas) thfile << theta << '\n';

        phfile << std::setprecision(15);
        for (const auto& phi : phis) phfile << phi << '\n';
    }

private:
    std::vector<interpPair> getInterpTheta(int, int);

    std::vector<interpPair> getInterpPhi(int, int);

    Map<vecXcd> getAlpha();

    HashMap<interpPair> getInterpPsi();

    void buildTranslationTable();

public:
    static std::vector<double> dists;
    static std::vector<vec3d> rhats;
    static std::array<vec3d, 316> dXs;

    // M2M interpolation tables
    std::vector<interpPair> interpTheta;
    std::vector<interpPair> interpPhi;

    // M2L translation table
    VecHashMap<arrXcd> transl;

    // Angular quantities
    std::vector<vec3d> khat;
    std::vector<mat23d> toThPh;
    std::vector<mat3d> ImRR;

    std::vector<double> thetas;  // theta samples
    std::vector<double> weights; // weights of theta samples
    std::vector<double> phis;    // phi samples
    int L;                       // M2L series truncation number

    int level; // level of this FMM level, with root at level 0
};
