#include "subrwg.h"

SubRWG::SubRWG(std::shared_ptr<Triangle> tri0, std::shared_ptr<Triangle> tri1)
    : RWG(nullptr, 0, std::move(tri0), std::move(tri1))
{
    using namespace Math;

    // If an RWG straddles two coarse mesh vertices, it does not contribute 
    // to the BC function at either vertex
    if (tris[0]->glIdxs[0] == tris[1]->glIdxs[0])
        glIdx = tris[0]->glIdxs[0];

    // Find common vertices
    int k = 0;
    for (int i = 0; i < 3; ++i) {
        const vec3d& X0 = tris[0]->Xs[i];
        for (int j = 0; j < 3; ++j) {
            const vec3d& X1 = tris[1]->Xs[j];
            if (vecEquals(X0,X1)) Xc[k++] = X0;
        }
    }
    const vec3d& dX = Xc[1]-Xc[0];
    center = (Xc[0]+Xc[1])/2.0;
    leng = dX.norm();

    /*
    // Find non-common vertices
    for (int i = 0; i < 2; ++i)
        for (const auto& X : tris[i]->Xs)
            if (!vecEquals(X,Xc[0]) && !vecEquals(X,Xc[1]))
                Xnc[i] = X;

    // Reorder tris if needed
    const vec3d& nhat0 = dX.cross(Xnc[0] - Xc[0]);
    const vec3d& nhat1 = dX.cross(Xnc[1] - Xc[0]);
    assert(nhat0.dot(nhat0 - nhat1) > 0);
    // std::cout << nhat0.dot(nhat0 - nhat1) << '\n';
    // if (nhat0.dot(nhat0 - nhat1) < 0) std::swap(tris[0],tris[1]);
    */
}

void SubRWG::buildVertsToSubrwgs(int numVerts) {
    // Remove duplicate subrwgs from global list based on their centers
    std::cout << glSubrwgs.size() << " subrwgs before removing duplicates\n";

    std::sort(glSubrwgs.begin(), glSubrwgs.end(), 
        [](std::shared_ptr<SubRWG> rwg0, std::shared_ptr<SubRWG> rwg1) {
            return Math::vecLessThan(rwg0->center,rwg1->center);
        }
    );

    glSubrwgs.erase(
        std::unique(glSubrwgs.begin(),glSubrwgs.end(),
        [](std::shared_ptr<SubRWG> rwg0,std::shared_ptr<SubRWG> rwg1) {
            return Math::vecEquals(rwg0->center,rwg1->center);
        }), 
        glSubrwgs.end()
    );

    std::cout << glSubrwgs.size() << " subrwgs after removing duplicates\n";

    // For each subrwg, map its coarse mesh vertex to itself
    vertsToSubrwgs.resize(numVerts);
    for (const auto& rwg : glSubrwgs) {
        if (rwg->glIdx)
            vertsToSubrwgs[rwg->glIdx.value()].push_back(std::move(rwg));
    }

    /*
    int vIdx = 0;
    for (const auto& rwgs : vertsToSubrwgs) {
        std::cout << vIdx++ << "\n";
        for (const auto& rwg : rwgs)
            std::cout << rwg->Xc[0] << ' ' << rwg->Xc[1] << '\n';
        std::cout << '\n';
    }
    */
}

void SubRWG::setOriented(const vec3d& Xref, const vec3d& nhat, const vec3d& ehat) {
    const auto& [X0,X1] = Xc;
    assert(Math::vecEquals(Xref,X0) || Math::vecEquals(Xref,X1));
    // std::cout << X_bc << ' ' << X0 << ' ' << X1 << '\n';
    const vec3d& rhat = (Math::vecEquals(Xref,X0) ? X1-X0 : X0-X1).normalized();
    const double angle = atan2(nhat.dot(ehat.cross(rhat)),ehat.dot(rhat));

    // std::cout << angle << '\n';

    oriented = (angle < 0.0 ? angle+2.0*PI : angle);
}