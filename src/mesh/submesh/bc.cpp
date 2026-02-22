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
        std::vector<SubRWG> subs;
        std::transform(iSubs.begin(), iSubs.end(), std::back_inserter(subs),
            [](int iSub) { return glSubrwgs[iSub]; } );

        // WLOG pick a surface normal (e.g. of first triangle of first subRWG)
        const vec3d& nhat = (subs[0].getTris())[0].nhat; 
        
        // Set CCW oriented angle for each subRWG relative to common edge
        for (auto& sub : subs) sub.setOriented(iVert,nhat,ehat);

        // Sort subRWGs by oriented angle
        std::sort(subs.begin(),subs.end(), 
            [](const SubRWG& rwg0, const SubRWG& rwg1) {
                return rwg0.oriented < rwg1.oriented;
            }
        );
        assert(approxZero(subs[0].oriented)); // first RWG should be along common edge

        // Recollect indices of sorted subRWGs at this vertex
        std::transform(subs.begin(), subs.end(), iSubs.begin(),
            [](const SubRWG& rwg) { return rwg.iSrc; } );

        // Discard RWG along common edge
        iSubs.erase(iSubs.begin()); 

        // Compute pcoeffs
        const int nSubs = iSubs.size();
        assert(nSubs%2); // for closed mesh, vertex should have odd # of subRWGs
        const double nTris = (nSubs+1.0)/2.0;

        std::vector<double> pcoeffs(iSubs.size());
        for (size_t j = 0; j < nSubs; ++j)
            pcoeffs[j] = (nTris - j - 1) / (2.0 * nTris);

        iSubss[idx] = std::move(iSubs);
        pcoeffss[idx++] = std::move(pcoeffs);
    }
};

// Receive rvals from fine RWGs and add to rval of this BC
void Mesh::BC::receiveRvals() {
    rval = 0.0;

    int iVert = 0;
    for (const auto& iSubs : iSubss) {
        for (int j = 0; j < iSubs.size(); ++j) {
            const auto& subrwg = glSubrwgs[iSubs[j]];

            rval += pcoeffss[iVert][j] // precompute fine RWG -> BC coeffs
                * subrwg.rval 
                * Math::sign(iVert)
                ;
        }

        ++iVert;
    }
}