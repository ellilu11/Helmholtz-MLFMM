#include "fmm.h"

void FMM::initGlobal(
    const Config& config_,
    const std::shared_ptr<Excitation::PlaneWave>& Einc,
    int nsrcs)
{
    config = config_;
    wavenum = Einc->wavenum;

    lvec = std::make_shared<vecXcd>(vecXcd::Zero(nsrcs));
    rvec = std::make_shared<vecXcd>(vecXcd::Zero(nsrcs));
    currents = std::make_shared<vecXcd>(vecXcd::Zero(nsrcs)); // assume I = 0 initially

    resetLeaves();
}

void FMM::buildTables() {
    // std::cout << "   (Lvl,Nth,Nph) =\n";
    angles.reserve(maxLevel+1);
    for (int level = 0; level <= maxLevel; ++level)
        angles.emplace_back(level);

    Tables::buildDists();
    tables.reserve(maxLevel+1);
    for (int level = 0; level <= maxLevel; ++level)
        tables.emplace_back(level, maxLevel);
}
