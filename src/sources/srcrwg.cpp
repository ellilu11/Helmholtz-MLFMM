#include "srcrwg.h"

SrcRWG::SrcRWG(
    std::shared_ptr<Excitation::PlaneWave> Einc,
    size_t iSrc,
    const Eigen::Vector4i& idx4)
    : RWG(std::move(Einc), iSrc, idx4)
{
    iCenter = Triangle::glEdgeToMid.at(
        makeUnordered(iVertsC[0], iVertsC[1]));

    buildVoltage(); 

    //std::cout << "Built srcRWG #" << iSrc << " w/ common vertices # " 
    //    << iVertsC[0] << ' '<< iVertsC[1] << " and non-common vertices # "
    //    << iVertsNC[0] << ' ' << iVertsNC[1] << " and center # "
    //    << iCenter << "\n";
};

void SrcRWG::buildSubIdx() {
    auto getIdxMid = [&](int i0, int i1) {
        return Triangle::glEdgeToMid.at(makeUnordered(i0, i1));
    };

    auto getIdxSub = [&](int i0, int i1) {
        return SubRWG::glEdgeToSub.at(makeUnordered(i0, i1));
    };

    int iPair = 0;
    for (auto iTri : iTris) {
        const int iCenter = Triangle::glTris[iTri].iCenter; // center of coarse tri
        const int iMid0 = getIdxMid(iVertsNC[iPair], iVertsC[0]); // midpoint between non-common and 0th common vertex
        const int iMid1 = getIdxMid(iVertsNC[iPair], iVertsC[1]); // midpoint between non-common and 1st common vertex

        const int iPair5 = 5*iPair;

        iSubs[iPair5] = getIdxSub(iMid0, iCenter); // edge 2 or 12
        iSubs[iPair5+1] = getIdxSub(iMid1, iCenter); // edge 3 or 14
        iSubs[iPair5+2] = getIdxSub(iVertsC[0], iCenter); // edge 4 or 9
        iSubs[iPair5+3] = getIdxSub(iVertsC[1], iCenter); // edge 6 or 11
        iSubs[iPair5+4] = getIdxSub(iVertsC[iPair], this->iCenter); // edge 7 or 8

        ++iPair;
    }

    //std::cout << "SrcRWG #" << iSrc << " contains subRWGs # ";
    //for (auto i : iSubs) std::cout << i << ' ';
    //std::cout << '\n';
}

/*
void SrcRWG::buildAncestry() {
    // Find all subedges of this RWG
    const auto& edgeMap = SubRWG::glEdgeToSub;

    for (auto i : iTris) {
        for (int j = 0; j < 6; ++j) {
            const int iSub = Triangle::NVerts + 6*i + j;
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
        const auto& iSubtris = Triangle::glEdgeToTri[edge];
        
        const int iBaseTri0 = (iSubtris[0] - Triangle::NTris) / 6;
        const int iBaseTri1 = (iSubtris[1] - Triangle::NTris) / 6;

        if (iBaseTri0 == iTris[0] || iBaseTri0 == iTris[1] || iBaseTri1 == iTris[0] || iBaseTri1 == iTris[1])
            iBases.push_back(iSrc);
    }//
}
*/

/*
void SrcRWG::buildSubRWGs() {
    const Eigen::MatrixXi subIdxs{
        {0, 0, 4, 1, 2, 3},
        {5, 1, 5, 2, 3, 4}
    };

    int iTri = 0;
    std::vector<Triangle> midSubtris;
    
    for (const auto& tri : tris) {
        const auto& subtris = tri->getSubtris(vec3i(iVertsNC[iTri],iVertsC[0],iVertsC[1]));
        const int iTri8 = 8*iTri;

        for (int iSub = 0; iSub < 6; ++iSub) {
            const auto& subtris0 = subtris[subIdxs(0,iSub)];
            const auto& subtris1 = subtris[subIdxs(1,iSub)];

            auto subrwg = std::make_shared<SubRWG>(subtris0,subtris1);

            subrwgs[iTri8+iSub] = subrwg; // TODO: Move assign
            glSubrwgs.push_back(std::move(subrwg));

            if (iSub == 2 || iSub == 3) midSubtris.push_back(subtris[iSub]);
        }

        ++iTri;
    }

    // Construct RWGs along common edge of parent RWG
    for (int iSub = 0; iSub < 2; ++iSub) {
        const auto& subtris0 = midSubtris[iSub];
        const auto& subtris1 = midSubtris[iSub+2];

        auto subrwg = std::make_shared<SubRWG>(subtris0,subtris1);
        subrwgs[iSub+6] = subrwg; // TODO: Move assign
        glSubrwgs.push_back(std::move(subrwg));
    }

    //for (const auto& rwg : subrwgs) {
    //    for (auto& idx : rwg->glIdxs)
    //        std::cout << idx << ' ';
    //    std::cout << '\n';
    //}
        // std::cout << rwg->Xc[0] << ' ' << rwg->Xc[1] << '\n';
}
*/