#include "bc.h"

Mesh::BC::BC(SrcRWG* const rwg) : base(rwg) {
    using namespace Math;

    const auto& vertsC = base->getVertsC();
    vec3d ehat = (vertsC[1]-vertsC[0]).normalized();

    int idx = 0;
    for (auto iVert : base->iVertsC) {
        const auto& vert = vertsC[idx];
        auto vertrwgs = vertsToSubrwgs[iVert];

        /* Find average nhat at X_bc
        vec3d nhat = vec3d::Zero();
        for (const auto& rwg : vertrwgs)
            for (const auto& tri : rwg->tris) {
                nhat += tri->alpha * tri->nhat;
                // std::cout << tri->alpha * tri->nhat << '\n';
            }
        nhat.normalize();
        //std::cout << nhat << '\n';
        */

        ehat *= sign(idx);
        const vec3d& nhat = (vertrwgs[0].getTris())[0].nhat; // WLOG pick an nhat
        
        for (auto& rwg : vertrwgs) {
            rwg.setOriented(vert,nhat,ehat);
            std::cout << rwg.oriented << ' ';
        }
        std::cout << '\n';

        // std::cout << "Vertex " << iVert << " has " << vertrwgs.size() << " subrwgs\n";

        std::sort(vertrwgs.begin(),vertrwgs.end(), 
            [](const SubRWG& rwg0, const SubRWG& rwg1) {
                return rwg0.oriented < rwg1.oriented;
            }
        );

        numRWGs[idx++] = vertrwgs.size();

        subrwgs.insert(subrwgs.end(),vertrwgs.begin(),vertrwgs.end());

        //for (const auto& rwg : vertrwgs) std::cout << rwg->oriented << ' ';
        // std::cout << '\n';
    }

    std::cout << '\n';
};