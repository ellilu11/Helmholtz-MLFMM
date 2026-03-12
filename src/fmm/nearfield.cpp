#include "nearfield.h"
#include "../mesh/tripair.h"

FMM::Nearfield::Nearfield(int nsrcs) 
    : nearMat(nsrcs, nsrcs)
{
    std::cout << " Building nearfield matrix...     ";

    auto start = Clock::now();

    std::vector<Eigen::Triplet<cmplx>> trips;

    findNodePairs();
    buildTriPairs();

    buildPairRads(trips);
    buildSelfRads(trips);

    nearMat.setFromTriplets(trips.begin(), trips.end());
    nearMat.makeCompressed();

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

void FMM::Nearfield::buildPairRads(std::vector<Eigen::Triplet<cmplx>>& trips) {
    // Precompute trips capacity
    size_t totalCap = 0;
    for (const auto& nearPair : nearPairs) {
        const auto [obsLeaf, srcNode] = nearPair.pair;
        size_t nObss = obsLeaf->srcs.size(), nSrcs = srcNode->srcs.size();
        totalCap += 2*nObss*nSrcs;
    }
    trips.reserve(trips.capacity() + totalCap);

    for (auto& nearPair : nearPairs) {
        const auto [obsLeaf, srcNode] = nearPair.pair;
        // assert(obsLeaf < srcNode);

        size_t nObss = obsLeaf->srcs.size(), nSrcs = srcNode->srcs.size();
        nearPair.rads.resize(nObss*nSrcs);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nObss; ++iObs) {
            auto obs = obsLeaf->srcs[iObs];
            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                auto src = srcNode->srcs[iSrc];

                cmplx rad = Phys::C * config.k * obs->getIntegratedRad(src);

                nearPair.rads[pairIdx++] = rad;
                trips.emplace_back(obs->getIdx(), src->getIdx(), rad);
                trips.emplace_back(src->getIdx(), obs->getIdx(), rad);
            }
        }
    }
}

void FMM::Nearfield::buildSelfRads(std::vector<Eigen::Triplet<cmplx>>& trips) {
    // Precompute trips capacity
    size_t totalCap = 0;
    for (const auto& selfPair : selfPairs) {
        const auto [leaf, srcLeaf] = selfPair.pair;
        size_t nSrcs = leaf->srcs.size();
        totalCap += nSrcs*nSrcs;
    }
    trips.reserve(trips.capacity() + totalCap);

    for (auto& selfPair : selfPairs) {
        const auto [leaf, srcLeaf] = selfPair.pair;
        assert(leaf == srcLeaf);

        size_t nSrcs = leaf->srcs.size();
        selfPair.rads.resize(nSrcs*(nSrcs+1)/2);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nSrcs; ++iObs) { // iObs = 0
            auto obs = leaf->srcs[iObs];

            for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) { // iSrc <= iObs 
                auto src = leaf->srcs[iSrc];

                cmplx rad = Phys::C * config.k * obs->getIntegratedRad(src);

                selfPair.rads[pairIdx++] = rad;
                trips.emplace_back(obs->getIdx(), src->getIdx(), rad);
                if (iSrc != iObs) // Only add self-term contribution once!
                    trips.emplace_back(src->getIdx(), obs->getIdx(), rad);
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

            cmplx rad = nearPair.rads[pairIdx++];

            rvecObs += lvecSrc * rad;
            rvecSrc += lvecObs * rad;
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

            cmplx rad = selfPair.rads[pairIdx++];

            // Only add self-term contribution once!
            if (iSrc == iObs) {
                rvecObs += lvecSrc * rad;
                continue;
            }

            rvecObs += lvecSrc * rad;
            rvecSrc += lvecObs * rad;
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
             cmplx rad = Phys::C * config.k * rads[iTri+iSrc];

             mat(iObs, iSrc) = rad;
             mat(iSrc, iObs) = rad;
         }
     }

     return mat;
 }

/*
sparseMat<cmplx> FMM::Nearfield::getNearMatrix(int nsrcs) const {


    // Near pairs first
    for (const auto& nearPair : nearPairs) {
        const auto& [obsLeaf, srcLeaf] = nearPair.pair;

        const SrcVec& srcs = obsLeaf->srcs;
        const SrcVec& srcSrcs = srcLeaf->srcs;
        size_t nObs = srcs.size(), nSrcs = srcSrcs.size();

        triplets.reserve(triplets.capacity() + 2*nObs*nSrcs); // symmetric entries

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nObs; ++iObs) {
            size_t obsIdx = srcs[iObs]->getIdx();

            for (size_t iSrc = 0; iSrc < nSrcs; ++iSrc) {
                size_t srcIdx = srcSrcs[iSrc]->getIdx();

                cmplx rad = Phys::C * config.k * nearPair.rads[pairIdx++];
                triplets.emplace_back(obsIdx, srcIdx, rad);
                triplets.emplace_back(srcIdx, obsIdx, rad);
            }
        }
     }

    // Then self pairs
     for (const auto& selfPair : selfPairs) {
        const auto& [leaf, srcLeaf] = selfPair.pair;
        assert(leaf == srcLeaf);
        
        const SrcVec& srcs = leaf->srcs;
        size_t nSrcs = srcs.size();

        triplets.reserve(triplets.capacity() + nSrcs*(nSrcs+1)/2);

        int pairIdx = 0;
        for (size_t iObs = 0; iObs < nSrcs; ++iObs) {
            size_t obsIdx = srcs[iObs]->getIdx();
            
            for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) {
                size_t srcIdx = srcs[iSrc]->getIdx();
                
                cmplx rad = Phys::C * config.k * selfPair.rads[pairIdx++];

                triplets.emplace_back(obsIdx, srcIdx, rad);
                if (iSrc != iObs)
                    triplets.emplace_back(srcIdx, obsIdx, rad);
            }
        }
     }

     mat.setFromTriplets(triplets.begin(), triplets.end());
     mat.makeCompressed();

     // Print mat for debugging
     std::ofstream ofs("out/nfmat.txt");
     matXcd denseMat(mat);
     for (int i = 0; i < denseMat.rows(); ++i) {
         for (int j = 0; j < denseMat.cols(); ++j) {
             ofs << denseMat(i, j).real() << " ";
         }
         ofs << "\n";
     }//

     return mat;
 }
 */
