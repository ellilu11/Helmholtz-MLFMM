#include "srcrwg.h"

Mesh::SrcRWG::SrcRWG(
    std::shared_ptr<Excitation::PlaneWave> Einc,
    size_t iSrc,
    const Eigen::Vector4i& idx4)
    : RWG(std::move(Einc), iSrc, idx4)
{
    buildVoltage(); 

    //std::cout << "Built srcRWG #" << iSrc << " w/ common vertices # " 
    //    << iVertsC[0] << ' '<< iVertsC[1] << " and non-common vertices # "
    //    << iVertsNC[0] << ' ' << iVertsNC[1] << "\n";
};

void Mesh::SrcRWG::buildSubIdx() {
    auto getIdxMid = [&](int i0, int i1) {
        return edgeToMid.at(makeUnordered(i0, i1));
    };

    auto getIdxSub = [&](int i0, int i1) {
        return fineEdgeToSub.at(makeUnordered(i0, i1));
    };

    auto iMid1 = getIdxMid(iVertsC[0], iVertsC[1]); // midpoint of common edge

    int iPair = 0;
    for (auto iTri : iTris) {
        auto iCenter = glTris[iTri].iCenter; // center of coarse tri
        auto iMid0 = getIdxMid(iVertsNC[iPair], iVertsC[0]); // midpoint of (non-common, 0th common)
        auto iMid2 = getIdxMid(iVertsNC[iPair], iVertsC[1]); // midpoint of (non-common, 1st common)

        auto i5 = 5*iPair;

        // Use a loop?
        iSubs[i5] = getIdxSub(iMid0, iCenter); // edge 2 or 12
        iSubs[i5+1] = getIdxSub(iMid2, iCenter); // edge 3 or 14
        iSubs[i5+2] = getIdxSub(iVertsC[0], iCenter); // edge 4 or 9
        iSubs[i5+3] = getIdxSub(iVertsC[1], iCenter); // edge 6 or 11
        iSubs[i5+4] = getIdxSub(iVertsC[iPair], iMid1); // edge 7 or 8

        ++iPair;
    }

    //std::cout << "SrcRWG #" << iSrc << " contains subRWGs # ";
    //for (auto i : iSubs) std::cout << i << ' ';
    //std::cout << '\n';

    /*
    for (auto i : iSubs) {
        //const auto& subrwg = SubRWG::glSubrwgs[i];
        //std::cout << subrwg.iVertsC[0] << ' ' << subrwg.iVertsC[1] << '\n';
         const auto& vertsC = glSubrwgs[i].getVertsC();
         std::cout << vertsC[0] << " " << vertsC[1] << '\n';
    }
    */
}

/*
void SrcRWG::buildAncestry() {
    // Find all subedges of this RWG
    const auto& edgeMap = SubRWG::fineEdgeToSub;

    for (auto i : iTris) {
        for (int j = 0; j < 6; ++j) {
            const int iSub = Triangle::nverts + 6*i + j;
            const auto& subtri = Triangle::glTris[iSub];

            for (int k = 0; k < 3; ++k) {
                const pair2i edge = makeUnordered(
                    subtri.iVerts[k],
                    subtri.iVerts[(k+1)%3]
                );

                // TODO: Exclude edges at boundary of this RWG
                const int iIdxSub = edgeMap.at(edge);
                iSubs.push_back(iIdxSub);
            }
        }
    }

    assert(iSubs.size() == 14);

    //
    for (const auto& edge : subedges) {
        const auto& iSubtris = Triangle::fineEdgeToTri[edge];
        
        const int iBaseTri0 = (iSubtris[0] - Triangle::ntris) / 6;
        const int iBaseTri1 = (iSubtris[1] - Triangle::ntris) / 6;

        if (iBaseTri0 == iTris[0] || iBaseTri0 == iTris[1] || iBaseTri1 == iTris[0] || iBaseTri1 == iTris[1])
            iBases.push_back(iSrc);
    }//
}
*/