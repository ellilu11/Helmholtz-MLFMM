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

void FMM::buildRadPats() {
    for (const auto& leaf : leaves)
        leaf->buildRadPats();
}

/*
void FMM::buildRadPats() {
    // Group leaves by level
    std::vector<NodeVec> leveledLeaves(maxLevel+1);
    for (const auto& leaf : leaves)
        leveledLeaves[leaf->getLevel()].push_back(leaf);
    
    for (int level = 0; level <= maxLevel; ++level) {
        const auto& leaves = leveledLeaves[level];
        if (leaves.empty()) continue;

        for (const auto& leaf : leveledLeaves[level])
            leaf->buildRadPats();
    }
}*/

void FMM::evaluateSols() {
    auto start = Clock::now();
    for (const auto& leaf : leaves)
        leaf->evalFarSols();
    t.L2T += Clock::now() - start;
}

void FMM::addInterpCoeffs(
    const Coeffs& inCoeffs, Coeffs& outCoeffs, int srcLvl, int tgtLvl)
{
    int order = config.interpOrder;

    auto [mth, mph] = angles[srcLvl].getNumAngles();
    auto [nth, nph] = angles[tgtLvl].getNumAngles();
    assert(!(mph%2)); // mph needs to be even

    int tblLvl = std::min(srcLvl, tgtLvl);
    const auto& interpTheta = tables[tblLvl].interpTheta;
    const auto& interpPhi = tables[tblLvl].interpPhi;

    // Interpolate over theta
    Coeffs innerCoeffs(nth*mph);
    for (int jth = 0; jth < nth; ++jth) {
        const auto [interp, nearIdx] = interpTheta[jth];

        for (int ith = nearIdx+1-order; ith <= nearIdx+order; ++ith) {
            int ith_flipped = Math::flipIdxToRange(ith, mth);
            int k = ith - (nearIdx+1-order);

            bool outOfRange = ith != ith_flipped; // jth < 0 || jth >= mth;

            for (int iph = 0; iph < mph; ++iph) {
                int iph_shifted = iph;
                if (outOfRange) iph_shifted += ((iph < mph/2) ? mph/2 : -mph/2);

                size_t idxInner = jth*mph+iph, idxIn = ith_flipped*mph+iph_shifted;

                innerCoeffs.theta[idxInner] +=
                    interp[k] * inCoeffs.theta[idxIn] * Math::sign(outOfRange);
                innerCoeffs.phi[idxInner] +=
                    interp[k] * inCoeffs.phi[idxIn] * Math::sign(outOfRange);
            }
        }
    }

    // Interpolate over phi
    for (int jph = 0; jph < nph; ++jph) {
        const auto [interp, nearIdx] = interpPhi[jph];

        for (int iph = nearIdx+1-order; iph <= nearIdx+order; ++iph) {
            int iph_wrapped = Math::wrapIdxToRange(iph, mph);
            int k = iph - (nearIdx+1-order);

            for (int jth = 0; jth < nth; ++jth) {
                size_t idxOut = jth*nph+jph, idxInner = jth*mph+iph_wrapped;

                outCoeffs.theta[idxOut] += interp[k] * innerCoeffs.theta[idxInner];
                outCoeffs.phi[idxOut] += interp[k] * innerCoeffs.phi[idxInner];
            }
        }
    }
}

void FMM::addAnterpCoeffs(
    const Coeffs& inCoeffs, Coeffs& outCoeffs, int srcLvl, int tgtLvl)
{
    int order = config.interpOrder;

    auto [mth, mph] = angles[srcLvl].getNumAngles();
    auto [nth, nph] = angles[tgtLvl].getNumAngles();
    assert(!(nph%2)); // nph needs to be even

    const int tblLvl = std::min(srcLvl, tgtLvl);
    const auto& interpTheta = tables[tblLvl].interpTheta;
    const auto& interpPhi = tables[tblLvl].interpPhi;

    // Anterpolate over extended phi
    Coeffs innerCoeffs(mth*nph);
    for (int iph = 0; iph < mph; ++iph) {
        const auto [interp, nearIdx] = interpPhi[iph];

        for (int jph = -order; jph < nph+order; ++jph) {
            int k = jph - (nearIdx+1-order);

            // If iph \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            int jph_wrapped = Math::wrapIdxToRange(jph, nph);

            for (int ith = 0; ith < mth; ++ith) {
                size_t idxInner = ith*nph+jph_wrapped, idxIn = ith*mph+iph;

                innerCoeffs.theta[idxInner] += interp[k] * inCoeffs.theta[idxIn];
                innerCoeffs.phi[idxInner] += interp[k] * inCoeffs.phi[idxIn];
            }
        }
    }

    // Anterpolate over extended theta
    for (int ith = 0; ith < mth; ++ith) {
        const auto [interp, nearIdx] = interpTheta[ith];

        for (int jth = -order; jth < nth+order; ++jth) {
            int k = jth - (nearIdx+1-order);

            // If ith \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            int jth_flipped = Math::flipIdxToRange(jth, nth);

            bool outOfRange = jth != jth_flipped; // jth < 0 || jth >= nth;

            for (int jph = 0; jph < nph; ++jph) {
                int jph_shifted = jph;
                if (outOfRange) jph_shifted += ((jph < nph/2) ? nph/2 : -nph/2);

                size_t idxOut = jth_flipped*nph+jph_shifted, idxInner = ith*nph+jph;

                outCoeffs.theta[idxOut] +=
                    interp[k] * innerCoeffs.theta[idxInner] * Math::sign(outOfRange);
                outCoeffs.phi[idxOut] +=
                    interp[k] * innerCoeffs.phi[idxInner] * Math::sign(outOfRange);
            }
        }
    }
}