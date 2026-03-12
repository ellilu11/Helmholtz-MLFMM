#include "nearfield.h"
#include "../mesh/tripair.h"

FMM::Nearfield::Nearfield(size_t nsrcs) 
    : nearMat(nsrcs, nsrcs)
{
    std::cout << " Building nearfield matrix...     ";
    auto start = Clock::now();

    findNodePairs();
    buildTriPairs();

    buildNearMatrix();
    Mesh::glTriPairs.clear();

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";
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
        const auto& [leaf, srcLeaf] = selfPair;
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
        const auto& [obsLeaf, srcNode] = nearPair;
        const auto &iTris0 = obsLeaf->iTris, &iTris1 = srcNode->iTris;

        for (auto iTri0 : iTris0)
            for (auto iTri1 : iTris1) {
                pair2i pair = std::minmax(iTri0, iTri1);
                Mesh::glTriPairs.emplace(pair, Mesh::TriPair(pair));
            }
    }
}

/* getNearCapacity()
 * Get number of nonzero entries in near matrix to reserve space for triplets
 */
size_t FMM::Nearfield::getNearCapacity() {
    size_t cap = 0;

    for (const auto& nearPair : nearPairs) {
        const auto [obsLeaf, srcNode] = nearPair;
        size_t nObss = obsLeaf->srcs.size(), nSrcs = srcNode->srcs.size();
        nearPair.rads.resize(nObss*nSrcs);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nObss; ++iObs) {
            auto obs = obsLeaf->srcs[iObs];
            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                auto src = srcNode->srcs[iSrc];

                nearPair.rads[pairIdx++] = obs->getIntegratedRad(src);
            }
        }
    }
}

    for (const auto& selfPair : selfPairs) {
        const auto [leaf, srcLeaf] = selfPair;
        size_t nSrcs = leaf->srcs.size();
        selfPair.rads.resize(nSrcs*(nSrcs+1)/2);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nSrcs; ++iObs) { // iObs = 0
            auto obs = leaf->srcs[iObs];

            for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
                auto src = leaf->srcs[iSrc];

                selfPair.rads[pairIdx++] = obs->getIntegratedRad(src);
            }
        }
    }
}

/* buildNearMatrix()
 * From list of node pairs, build near matrix by computing pairwise contributions
 * between sources in each node pair and adding to nearMat
 */
void FMM::Nearfield::buildNearMatrix() {
    std::vector<Eigen::Triplet<cmplx>> trips;
    trips.reserve(getNearCapacity());

    // Build pair-node contributions to near matrix
    for (auto& nearPair : nearPairs) {
        const auto [obsLeaf, srcNode] = nearPair;
        size_t nObss = obsLeaf->srcs.size(), nSrcs = srcNode->srcs.size();

        for (size_t iObs = 0; iObs < nObss; ++iObs) {
            auto obs = obsLeaf->srcs[iObs];
            size_t obsIdx = obs->getIdx(); // global index of obs

            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                auto src = srcNode->srcs[iSrc];
                size_t srcIdx = src->getIdx(); // global index of src

            cmplx rad = Phys::C * config.k * nearPair.rads[pairIdx++];

            rvecObs += lvecSrc * rad;
            rvecSrc += lvecObs * rad;
        }
    }
}

    // Build self-node contributions to near matrix
    for (auto& selfPair : selfPairs) {
        const auto [leaf, srcLeaf] = selfPair;
        assert(leaf == srcLeaf);

        size_t nSrcs = leaf->srcs.size();

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nSrcs; ++iObs) { // iObs = 0
            auto obs = leaf->srcs[iObs];
            size_t obsIdx = obs->getIdx(); // global index of obs

            for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
                auto src = leaf->srcs[iSrc];
                size_t srcIdx = src->getIdx(); // global index of src

            cmplx rad = Phys::C * config.k * selfPair.rads[pairIdx++];

            // Only add self-term contribution once!
            if (iSrc == iObs) {
                rvecObs += lvecSrc * rad;
                continue;
            }

            rvecObs += lvecSrc * rad;
            rvecSrc += lvecObs * rad;
        }
    }

    nearMat.setFromTriplets(trips.begin(), trips.end());
    nearMat.makeCompressed();
}

/* evaluateSols()
 * Multiply near matrix by lvec and add to rvec to get nearfield contribution to rvec
 */ 
void FMM::Nearfield::evaluateSols() {
    auto start = Clock::now();
  
    assert(nearMat.cols() == Solver::lvec.rows());
    assert(nearMat.rows() == Solver::rvec.rows());

    Solver::rvec += nearMat * Solver::lvec;

    t.S2T += Clock::now() - start;
}

/* printNearMatrix()
 * Print near matrix to file for debugging
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
             cmplx rad = Phys::C * config.k * rads[iTri+iSrc];

             mat(iObs, iSrc) = rad;
             mat(iSrc, iObs) = rad;
         }
     }

     return mat;
 }