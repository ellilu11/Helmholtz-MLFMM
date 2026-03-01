#include "nearfield.h"
#include "../mesh/tripair.h"

FMM::Nearfield::Nearfield() {
    findNodePairs();

    auto start = Clock::now();
    buildTriPairs();
    Time duration_ms = Clock::now() - start;
    std::cout << "  Built triangle pairs in " << duration_ms.count() << " ms\n";
    std::cout << "  # tri pairs: " << Mesh::glTriPairs.size() << '\n';

    //for (const auto& triPair : Mesh::glTriPairs)
    //    std::cout << "("
    //        << triPair.first.first << ","
    //        << triPair.first.second << ") ";
    //std::cout << '\n';

    start = Clock::now();
    buildPairRads();
    buildSelfRads();
    duration_ms = Clock::now() - start;
    std::cout << "  Built nearfield rads in " << duration_ms.count() << " ms\n";

    Mesh::glTriPairs.clear();
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
     
        const auto 
            &iTris0 = obsLeaf->iTris,
            &iTris1 = srcNode->iTris;

        for (auto iTri0 : iTris0)
            for (auto iTri1 : iTris1) {
                pair2i pair(iTri0, iTri1);
                Mesh::glTriPairs.emplace(pair, Mesh::TriPair(pair));
            }
    }
}

void FMM::Nearfield::buildPairRads() {
    for (auto& nearPair : nearPairs) {
        const auto [obsLeaf, srcNode] = nearPair.pair;
        // assert(obsLeaf < srcNode);

        size_t nObss = obsLeaf->srcs.size(), nSrcs = srcNode->srcs.size();
        nearPair.rads.resize(nObss*nSrcs);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nObss; ++iObs) {
            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                const auto obs = obsLeaf->srcs[iObs], src = srcNode->srcs[iSrc];

                nearPair.rads[pairIdx++] = obs->getIntegratedRad(src);
            }
        }
    }
}

void FMM::Nearfield::buildSelfRads() {
    for (auto& selfPair : selfPairs) {
        const auto [leaf, srcLeaf] = selfPair.pair;
        assert(leaf == srcLeaf);

        size_t nSrcs = leaf->srcs.size();
        selfPair.rads.resize(nSrcs*(nSrcs+1)/2);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < leaf->srcs.size(); ++iObs) { // iObs = 0
            const auto& obs = leaf->srcs[iObs];

            for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
                const auto& src = leaf->srcs[iSrc];

                selfPair.rads[pairIdx++] = obs->getIntegratedRad(src);
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

    const auto& srcs = obsLeaf->srcs;
    const auto& srcSrcs = srcLeaf->srcs;

    size_t nObs = srcs.size(), nSrcs = srcSrcs.size();

    std::vector<cmplx> solAtObss(nObs, 0.0);
    std::vector<cmplx> solAtSrcs(nSrcs, 0.0);

    int pairIdx = 0;
    for (size_t iObs = 0; iObs < nObs; ++iObs) {
        for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
            const auto obs = srcs[iObs], src = srcSrcs[iSrc];

            cmplx rad = nearPair.rads[pairIdx++];

            solAtObss[iObs] += states.lvec[src->getIdx()] * rad;
            solAtSrcs[iSrc] += states.lvec[obs->getIdx()] * rad;
        }
    }

    for (int n = 0; n < nObs; ++n)
        states.rvec[srcs[n]->getIdx()] += Phys::C * k * solAtObss[n];

    for (int n = 0; n < nSrcs; ++n)
        states.rvec[srcSrcs[n]->getIdx()] += Phys::C * k * solAtSrcs[n];
}

/* evalSelfSols()
 * (S2T) Evaluate sols at sources in this node due to other sources
 * in this node
 */
void FMM::Nearfield::evalSelfSols(const NearPair& selfPair) {
    const auto& [leaf, srcLeaf] = selfPair.pair;
    assert(leaf == srcLeaf);

    const auto& srcs = leaf->srcs;
    size_t nSrcs = srcs.size();

    std::vector<cmplx> solAtObss(nSrcs, 0.0);

    int pairIdx = 0;
    for (size_t iObs = 0; iObs < nSrcs; ++iObs) { // iObs = 0
        for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
            auto obs = srcs[iObs], src = srcs[iSrc];

            cmplx rad = selfPair.rads[pairIdx++];

            solAtObss[iObs] += states.lvec[src->getIdx()] * rad;
            solAtObss[iSrc] += states.lvec[obs->getIdx()] * rad;
        }
    }

    for (int n = 0; n < nSrcs; ++n)
        states.rvec[srcs[n]->getIdx()] += Phys::C * k * solAtObss[n];
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