#include "nearfield.h"
#include "../mesh/tripair.h"

FMM::Nearfield::Nearfield() {
    findNodePairs();
    buildTriPairs();

    buildPairRads();
    buildSelfRads();

    // Mesh::glTriPairs.clear();
}

/* findNodePairs()
 * From list of leaves, find all near neighbor leaf pairs
 */
void FMM::Nearfield::findNodePairs() {
    for (const auto& leaf : leaves) {
        selfPairs.emplace_back(leaf, leaf);

        for (const auto& nbor : leaf->nearNbors) {
            assert(nbor->isLeaf());
            if (leaf < nbor) nearPairs.emplace_back(leaf, nbor);
        }
    }

    for (const auto& pair : nonNearPairs) {
        const auto& [obsNode, srcNode] = pair;
        nearPairs.emplace_back(obsNode, srcNode);
    }
}

/* findTriPairs()
 * From list of node pairs, find all near neighbor triangle pairs
 * and populate Mesh::glTriPairs
 */
void FMM::Nearfield::buildTriPairs() {
    for (const auto& selfPair : selfPairs) {
        const auto& [leaf, srcLeaf] = selfPair.pair;
        assert(leaf == srcLeaf);
        const auto& iTris = leaf->iTris;

        for (auto iTri0 : iTris)
            for (auto iTri1 : iTris) {
                if (iTri0 > iTri1) continue;

                pair2i pair(iTri0, iTri1);
                Mesh::glTriPairs.emplace(pair, Mesh::TriPair(pair));
            }
    }

    for (const auto& nearPair : nearPairs) {
        const auto& [obsLeaf, srcNode] = nearPair.pair;
        const auto &iTris0 = obsLeaf->iTris, &iTris1 = srcNode->iTris;

        for (auto iTri0 : iTris0)
            for (auto iTri1 : iTris1) {
                pair2i pair = std::minmax(iTri0, iTri1);
                Mesh::glTriPairs.emplace(pair, Mesh::TriPair(pair));
            }
    }
}

void FMM::Nearfield::buildPairRads() {
    for (auto& nearPair : nearPairs) {
        const auto [obsLeaf, srcNode] = nearPair.pair;
        // assert(obsLeaf < srcNode);

        size_t nObss = obsLeaf->srcs.size(), nSrcs = srcNode->srcs.size();
        nearPair.efie.resize(nObss*nSrcs);
        nearPair.mfie.resize(nObss*nSrcs);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nObss; ++iObs) {
            auto obs = obsLeaf->srcs[iObs];
            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                auto src = srcNode->srcs[iSrc];

                nearPair.efie[pairIdx] = obs->getIntegratedEFIE(src);

                double mass = obs->getIntegratedMass(src);
                nearPair.mfie[pairIdx] =
                    { obs->getIntegratedMFIE(src)+mass, src->getIntegratedMFIE(obs)+mass };

                ++pairIdx;
            }
        }
    }

    std::cout << std::setprecision(3);
}

void FMM::Nearfield::buildSelfRads() {
    for (auto& selfPair : selfPairs) {
        const auto [leaf, srcLeaf] = selfPair.pair;
        assert(leaf == srcLeaf);

        size_t nSrcs = leaf->srcs.size();
        selfPair.efie.resize(nSrcs*(nSrcs+1)/2);
        selfPair.mfie.resize(nSrcs*(nSrcs+1)/2);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nSrcs; ++iObs) { // iObs = 0
            auto obs = leaf->srcs[iObs];

            for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
                auto src = leaf->srcs[iSrc];

                selfPair.efie[pairIdx] = obs->getIntegratedEFIE(src);

                double mass = obs->getIntegratedMass(src);
                selfPair.mfie[pairIdx] =
                    { obs->getIntegratedMFIE(src)+mass, src->getIntegratedMFIE(obs)+mass };

                ++pairIdx;
            }
        }
    }
}

/* (S2T) Evaluate sols at sources in this leaf due to sources in srcNode
 * and vice versa
 * srcNode : source node
 */
void FMM::Nearfield::evalPairSols(const NearPair& nearPair) {
    const auto& [obsLeaf, srcLeaf] = nearPair.pair;

    const SrcVec& srcs = obsLeaf->srcs;
    const SrcVec& srcSrcs = srcLeaf->srcs;
    size_t nObs = srcs.size(), nSrcs = srcSrcs.size();

    int pairIdx = 0;
    for (size_t iObs = 0; iObs < nObs; ++iObs) {
        auto obs = srcs[iObs];
        cmplx lvecObs = Solver::lvec[obs->getIdx()];
        cmplx& rvecObs = Solver::rvec[obs->getIdx()]; // &

        for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
            auto src = srcSrcs[iSrc];
            cmplx lvecSrc = Solver::lvec[src->getIdx()];
            cmplx& rvecSrc = Solver::rvec[src->getIdx()]; // &

            cmplx radAtObs = Phys::C * config.k *
                (config.alpha * nearPair.efie[pairIdx]
                    + config.beta * nearPair.mfie[pairIdx].first);
            cmplx radAtSrc = Phys::C * config.k *
                (config.alpha * nearPair.efie[pairIdx]
                    + config.beta * nearPair.mfie[pairIdx].second);
            ++pairIdx;

            rvecObs += lvecSrc * radAtObs;
            rvecSrc += lvecObs * radAtSrc;
        }
    }
}

/* evalSelfSols()
 * (S2T) Evaluate sols at sources in this node due to other sources
 * in this node
 */
void FMM::Nearfield::evalSelfSols(const NearPair& selfPair) {
    const auto& [leaf, srcLeaf] = selfPair.pair;
    assert(leaf == srcLeaf);

    const SrcVec& srcs = leaf->srcs;
    size_t nSrcs = srcs.size();

    int pairIdx = 0;
    for (size_t iObs = 0; iObs < nSrcs; ++iObs) { // iObs = 0
        auto obs = srcs[iObs];
        cmplx lvecObs = Solver::lvec[obs->getIdx()];
        cmplx& rvecObs = Solver::rvec[obs->getIdx()]; // &

        for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
            auto src = srcs[iSrc];
            cmplx lvecSrc = Solver::lvec[src->getIdx()];
            cmplx& rvecSrc = Solver::rvec[src->getIdx()]; // &

            cmplx radAtObs = Phys::C * config.k *
                (config.alpha * selfPair.efie[pairIdx]
                    + config.beta * selfPair.mfie[pairIdx].first);
            cmplx radAtSrc = Phys::C * config.k *
                (config.alpha * selfPair.efie[pairIdx]
                    + config.beta * selfPair.mfie[pairIdx].second);
            ++pairIdx;

            // Only add self-term contribution once!
            if (iSrc == iObs) {
                rvecObs += lvecSrc * radAtObs;
                continue;
            }

            rvecObs += lvecSrc * radAtObs;
            rvecSrc += lvecObs * radAtSrc;
        }
    }
}

/* evaluateSols()
 * Sum solutions at all sources in all leaves 
 */ 
void FMM::Nearfield::evaluateSols() {
    auto start = Clock::now();
  
    for (const auto& nearPair : nearPairs)
        evalPairSols(nearPair);

    for (const auto& selfPair : selfPairs)
        evalSelfSols(selfPair);

    t.S2T += Clock::now() - start;
}

/* getNearMatrix()
 * Get nearfield (impedance) matrix of all interactions between sources in this self pair
 */
 matXcd FMM::NearPair::getNearMatrix() const {
     const auto& [leaf, srcLeaf] = pair;
     assert(leaf == srcLeaf);

     const SrcVec& srcs = leaf->getSrcs();
     size_t nSrcs = srcs.size();
     matXcd mat(nSrcs, nSrcs);

     for (size_t iObs = 0; iObs < nSrcs; ++iObs) {
         size_t iTri = iObs*(iObs+1)/2; // triangular numbers

         for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) {
             cmplx radAtObs = Phys::C * config.k *
                 (config.alpha * efie[iTri+iSrc]
                     + config.beta * mfie[iTri+iSrc].first);
             cmplx radAtSrc = Phys::C * config.k *
                 (config.alpha * efie[iTri+iSrc]
                     + config.beta * mfie[iTri+iSrc].second);

             mat(iObs, iSrc) = radAtObs;
             mat(iSrc, iObs) = radAtSrc;
         }
     }

     return mat;
 }