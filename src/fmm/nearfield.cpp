#include "nearfield.h"
#include "../mesh/tripair.h"

FMM::Nearfield::Nearfield(size_t nsrcs) 
    : nearMat(nsrcs, nsrcs)
{
    std::cout << " Building nearfield matrix...     ";
    auto start = Clock::now();

    findNodePairs();
    buildTriPairs();

    std::cout << " # tripairs: " << Mesh::glTriPairs.size() << '\n';

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
    return
        std::accumulate(nearPairs.begin(), nearPairs.end(), 0,
            [](size_t sum, const auto& nearPair) {
                const auto& [obsLeaf, srcNode] = nearPair;
                size_t nObss = obsLeaf->srcs.size(), nSrcs = srcNode->srcs.size();
                return sum + 2*nObss*nSrcs;
            }
        ) +
        std::accumulate(selfPairs.begin(), selfPairs.end(), 0,
            [](size_t sum, const auto& selfPair) {
                const auto& [leaf, srcLeaf] = selfPair;
                size_t nSrcs = leaf->srcs.size();
                return sum + nSrcs*nSrcs;
            }
        );
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

                cmplx rad = Phys::C * config.k * obs->getIntegratedRad(src);

                trips.emplace_back(obsIdx, srcIdx, rad);
                trips.emplace_back(srcIdx, obsIdx, rad);
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

                cmplx rad = Phys::C * config.k * obs->getIntegratedRad(src);

                trips.emplace_back(obsIdx, srcIdx, rad);
                if (iSrc != iObs) // Only add self-term contribution once!
                    trips.emplace_back(srcIdx, obsIdx, rad);
            }
        }
    }

    nearMat.setFromTriplets(trips.begin(), trips.end());
    nearMat.makeCompressed();
}

/* evaluateSols()
 * (S2T) Multiply near matrix by lvec and add to rvec to get nearfield contribution to rvec
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
void FMM::Nearfield::printNearMatrix(const std::string& fname) const {
    std::ofstream of(fname);
    matXcd denseMat(nearMat);
    for (int i = 0; i < denseMat.rows(); ++i) {
        for (int j = 0; j < denseMat.cols(); ++j)
            of << denseMat(i, j).real() << ' ';
        of << '\n';
    }
}