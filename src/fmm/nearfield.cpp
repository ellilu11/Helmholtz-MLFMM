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

    // std::cout << "   # non-near pairs: " << nonNearPairs.size() << '\n';
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
        const auto& iTris0 = obsLeaf->iTris, & iTris1 = srcNode->iTris;

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

    // Build self-node contributions to near matrix
    for (auto& selfPair : selfPairs) {
        const auto [leaf, srcLeaf] = selfPair;
        assert(leaf == srcLeaf);

        size_t nSrcs = leaf->srcs.size();

        for (size_t iObs = 0; iObs < nSrcs; ++iObs) { // iObs = 0
            auto obs = leaf->srcs[iObs];
            size_t obsIdx = obs->getIdx(); // global index of obs

            for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
                auto src = leaf->srcs[iSrc];
                size_t srcIdx = src->getIdx(); // global index of src

                double mass = obs->getIntegratedMass(src);
                cmplx efie = config.C_efie * obs->getIntegratedEFIE(src);
                cmplx mfieObs = config.C_mfie * (obs->getIntegratedMFIE(src) + mass);

                assert(obsIdx < nearMat.rows() && srcIdx < nearMat.cols());
                trips.emplace_back(obsIdx, srcIdx, efie+mfieObs);

                if (iSrc != iObs) { // Only add self-term contribution once!
                    cmplx mfieSrc = config.C_mfie * (src->getIntegratedMFIE(obs) + mass);

                    assert(srcIdx < nearMat.rows() && obsIdx < nearMat.cols());
                    trips.emplace_back(srcIdx, obsIdx, efie+mfieSrc);
                }
            }
        }
    }

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

                double mass = obs->getIntegratedMass(src);
                cmplx efie = config.C_efie * obs->getIntegratedEFIE(src),
                    mfieObs = config.C_mfie * (obs->getIntegratedMFIE(src) + mass),
                    mfieSrc = config.C_mfie * (src->getIntegratedMFIE(obs) + mass);

                assert(obsIdx < nearMat.rows() && srcIdx < nearMat.cols());
                trips.emplace_back(obsIdx, srcIdx, efie+mfieObs);
                assert(srcIdx < nearMat.rows() && obsIdx < nearMat.cols());
                trips.emplace_back(srcIdx, obsIdx, efie+mfieSrc);
            }
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
void FMM::Nearfield::printNearMatrix(const std::string& fname) const {
    std::ofstream of(fname);
    matXcd denseMat(nearMat);
    for (int i = 0; i < denseMat.rows(); ++i) {
        for (int j = 0; j < denseMat.cols(); ++j)
            of << denseMat(i, j).real() << ' ';
        of << '\n';
    }
}


/* void FMM::Nearfield::evalPairSols(const NearPair& nearPair) {
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

            rvecObs += lvecSrc * nearPair.cfie[pairIdx].first;
            rvecSrc += lvecObs * nearPair.cfie[pairIdx].second;

            ++pairIdx;
        }
    }
}

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

            cmplx radAtObs = selfPair.cfie[pairIdx].first;
            cmplx radAtSrc = selfPair.cfie[pairIdx].second;
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
*/