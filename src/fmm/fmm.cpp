#include "fmm.h"

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
