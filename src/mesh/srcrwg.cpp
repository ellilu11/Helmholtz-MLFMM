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

void Mesh::SrcRWG::findSubRWGs() {
    auto getMidIdx = [&](int i0, int i1) {
        return edgeToMid.at(makeUnordered(i0, i1));
    };

    auto getSubIdx = [&](int i0, int i1) {
        return fineEdgeToSub.at(makeUnordered(i0, i1));
    };

    auto iMid1 = getMidIdx(iVertsC[0], iVertsC[1]); // midpoint of common edge

    int iPair = 0;
    for (auto iTri : iTris) {
        auto iCenter = glTris[iTri].iCenter; // center of coarse tri
        auto iMid0 = getMidIdx(iVertsNC[iPair], iVertsC[0]); // midpoint of (non-common, 0th common)
        auto iMid2 = getMidIdx(iVertsNC[iPair], iVertsC[1]); // midpoint of (non-common, 1st common)

        auto i5 = 5*iPair;

        // Use a loop?
        iSubs[i5] = getSubIdx(iMid0, iCenter);          // edge 2 or 12
        iSubs[i5+1] = getSubIdx(iMid2, iCenter);        // edge 3 or 14
        iSubs[i5+2] = getSubIdx(iVertsC[0], iCenter);   // edge 4 or 9
        iSubs[i5+3] = getSubIdx(iVertsC[1], iCenter);   // edge 6 or 11
        iSubs[i5+4] = getSubIdx(iVertsC[iPair], iMid1); // edge 7 or 8

        ++iPair;
    }

    /*
    for (auto i : iSubs) {
        //const auto& subrwg = SubRWG::glSubrwgs[i];
        //std::cout << subrwg.iVertsC[0] << ' ' << subrwg.iVertsC[1] << '\n';
         const auto& vertsC = glSubrwgs[i].getVertsC();
         std::cout << vertsC[0] << " " << vertsC[1] << '\n';
    }
    */
}

// Propagate rval to rval of fine RWGs
// TODO: Zero out fine RWG rvals in caller before calling this
void Mesh::SrcRWG::propagateRvals() {
    rval = states.rvec[iSrc]; // get updated rval from FMM

    for (size_t i = 0; i < iSubs.size(); ++i) {
        auto& subrwg = glSubrwgs[iSubs[i]];

        subrwg.addToRval(leng / subrwg.getLeng() * rcoeffs[i%5] * rval);
    }
}