#include "subrwg.h"

void Mesh::SubRWG::buildSubRWGs() {
    int iSub = 0;
    for (const auto& [edge, iTris] : glEdgeToTri) {
        const vec4i& idx4 = 
            { edge.first, edge.second, iTris[0], iTris[1]};

        if (iTris[1] < 0) continue; // Ignore edges not shared by two tris
        //std::cout << "Building SubRWG for edge ("
        //    << edge.first << ',' << edge.second << ") with tris # "
        //    << iTri[0] << ' ' << iTri[1] << '\n';

        glSubrwgs.emplace_back(iSub, idx4);
        glEdgeToSub.emplace(edge, iSub);
        ++iSub;
    }

    //for (const auto& [edge, iSub] : glEdgeToSub) {
    //    std::cout << "Edge (" << edge.first << ',' << edge.second
    //        << ") maps to subRWG # " << iSub << '\n';
    //}
}

Mesh::SubRWG::SubRWG(int iSub, const vec4i& idx4)
    : RWG(nullptr, iSub, idx4)
{
    const auto& [tri0, tri1] = getTris();

    // If an RWG straddles two coarse mesh vertices, 
    // it does not contribute to the BC function at either vertex
    if (tri0.iVerts[0] == tri1.iVerts[0])
        iVertsCoarse = tri0.iVerts[0];

    std::cout << "Built subRWG #" << iSrc << " w/ common vertices # "
        << iVertsC[0] << ' '<< iVertsC[1] << " and non-common vertices # "
        << iVertsNC[0] << ' ' << iVertsNC[1] << " and coarse vertex # "
        << (iVertsCoarse.has_value() ? std::to_string(iVertsCoarse.value()) : "N/A") << "\n";

    /* Reorder tris if needed
    const vec3d& nhat0 = dX.cross(Xnc[0] - Xc[0]);
    const vec3d& nhat1 = dX.cross(Xnc[1] - Xc[0]);
    assert(nhat0.dot(nhat0 - nhat1) > 0);
    // std::cout << nhat0.dot(nhat0 - nhat1) << '\n';
    // if (nhat0.dot(nhat0 - nhat1) < 0) std::swap(tris[0],tris[1]);
    */
}

void Mesh::SubRWG::buildVertsToSubrwgs(int numVerts) {
    // For each subrwg, map its coarse mesh vertex to itself
    vertsToSubrwgs.resize(numVerts);
    for (const auto& rwg : glSubrwgs) {
        if (rwg.iVertsCoarse)
            vertsToSubrwgs[rwg.iVertsCoarse.value()].push_back(std::move(rwg));
    }

    /*
    int vIdx = 0;
    for (const auto& vertToRwg : vertsToSubrwgs) {
        std::cout << "Vertex # " << vIdx++ << " has subRWGs with common vertices # ";
        for (const auto& rwg : vertToRwg)
            std::cout << '(' << rwg->iVertsC[0] << ',' << rwg->iVertsC[1] << ") ";
        std::cout << '\n';
    }
    */
}

void Mesh::SubRWG::setOriented(const vec3d& Xref, const vec3d& nhat, const vec3d& ehat) {
    const auto& [X0,X1] = getVertsC();
    assert(Math::vecEquals(Xref,X0) || Math::vecEquals(Xref,X1));
    // std::cout << X_bc << ' ' << X0 << ' ' << X1 << '\n';
    const vec3d& rhat = (Math::vecEquals(Xref,X0) ? X1-X0 : X0-X1).normalized();
    const double angle = atan2(nhat.dot(ehat.cross(rhat)),ehat.dot(rhat));

    // std::cout << angle << '\n';

    oriented = (angle < 0.0 ? angle+2.0*PI : angle);
}