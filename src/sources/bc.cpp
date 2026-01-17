#include "bc.h"

BC::BC(SrcRWG* const rwg) : base(rwg) {
    using namespace Math;

    const auto& Xc = base->Xc;

    int idx = 0;
    for (auto glIdx : base->glIdxc) {
        const auto& X_bc = Xc[idx];
        SubRWGVec vertrwgs = SubRWG::vertsToSubrwgs[glIdx]; // Find all subrwgs with X_bc as vertex

        /* Find average nhat at X_bc
        vec3d nhat = vec3d::Zero();
        for (const auto& rwg : vertrwgs)
            for (const auto& tri : rwg->tris) {
                nhat += tri->alpha * tri->nhat;
                // std::cout << tri->alpha * tri->nhat << '\n';
            }
        nhat.normalize();
        */
        //std::cout << nhat << '\n';

        const vec3d& ehat = (Xc[1]-Xc[0]).normalized() * sign(idx);
        const vec3d& nhat = vertrwgs[0]->tris[0]->nhat; // WLOG pick an nhat
        
        for (const auto& rwg : vertrwgs) {
            rwg->setOriented(X_bc,nhat,ehat);
            std::cout << rwg->oriented << ' ';
        }
        std::cout << '\n';

        std::cout << "Vertex " << glIdx << " has " << vertrwgs.size() << " subrwgs\n";

        std::sort(vertrwgs.begin(),vertrwgs.end(), 
            [](std::shared_ptr<SubRWG> rwg0, std::shared_ptr<SubRWG> rwg1) {
                return rwg0->oriented < rwg1->oriented;
            }
        );

        numRWGs[idx++] = vertrwgs.size();

        subrwgs.insert(subrwgs.end(),vertrwgs.begin(),vertrwgs.end());

        //for (const auto& rwg : vertrwgs) std::cout << rwg->oriented << ' ';
        // std::cout << '\n';
    }

    std::cout << '\n';
};