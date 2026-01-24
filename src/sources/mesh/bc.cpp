#include "bc.h"

Mesh::BC::BC(SrcRWG* const rwg) : base(rwg) {
    using namespace Math;

    const auto& vertsC = base->getVertsC();
    vec3d ehat = (vertsC[1]-vertsC[0]).normalized(); // common edge unit vector

    int idx = 0;
    for (auto iVert : base->iVertsC) {
        ehat *= sign(idx);
        const auto& vert = glVerts[iVert]; // vertsC[idx];
        auto iSubs = vertToSubs[iVert];

        // Collect subRWGs at this vertex
        std::vector<SubRWG> vertrwgs;
        std::transform(iSubs.begin(), iSubs.end(), std::back_inserter(vertrwgs),
            [](int iSub) { return glSubrwgs[iSub]; } );

        // WLOG pick a surface normal (e.g. of first triangle of first subRWG)
        const vec3d& nhat = (vertrwgs[0].getTris())[0].nhat; 
        
        // Set oriented angle for each subRWG relative to common edge
        for (auto& rwg : vertrwgs) rwg.setOriented(iVert,nhat,ehat);

        // Sort subRWGs by oriented angle
        std::sort(vertrwgs.begin(),vertrwgs.end(), 
            [](const SubRWG& rwg0, const SubRWG& rwg1) {
                return rwg0.oriented < rwg1.oriented;
            }
        );
        assert(approxZero(vertrwgs[0].oriented)); // 0th RWG should be along common edge

        // Recollect indices of sorted subRWGs at this vertex
        std::transform(vertrwgs.begin(), vertrwgs.end(), iSubs.begin(),
            [](const SubRWG& rwg) { return rwg.iSrc; } );

        // Discard RWG along common edge
        iSubs.erase(iSubs.begin()); 

        iSubss[idx++] = std::move(iSubs);
    }
};

// Receive rvals from fine RWGs and add to rval of this BC
void Mesh::BC::receiveRvals() {
    rval = 0.0;

    int iVert = 0;
    for (const auto& iSubs : iSubss) {
        assert(iSubs.size()%2); // for closed mesh, vertex should have odd number of subRWGs
        const double ntris = (iSubs.size()+1.0)/2.0;

        for (int j = 0; j < iSubs.size(); ++j) {
            const auto& subrwg = glSubrwgs[iSubs[j]];

            rval += (ntris-j-1) / (2.0*ntris)
                * subrwg.rval 
                * Math::sign(iVert)
                ;
        }

        ++iVert;
    }
}