#include "bc.h"

BC::BC(RWG* rwg) :
    X0(rwg->X0), X1(rwg->X1)
{
    /*
    const vec3d& ehat = (X1-X0).normalized();

    // Find all subtris with X0 as vertex
    TriVec subtris0;
    const int idx0 = rwg->idx0;
    for (auto iTri : Triangle::vertsToSubtris[idx0])
        subtris0.push_back(Triangle::glSubtris[iTri]);

    std::cout << "BC has " << subtris0.size() << " subtris centered at vertex #" << idx0 << "\n";

    // Sort subtris by oriented angle relative to RWG common edge
    auto oriented = [](const vec3d& nhat, const vec3d& ehat, const vec3d& rhat) {
        double angle = atan2(nhat.dot(ehat.cross(rhat)), ehat.dot(rhat));
        // double angle = atan2((ehat.cross(rhat)).norm(), ehat.dot(rhat));
        return (angle < 0.0 ? angle+2.0*PI : angle);
    };


    vec3d nhatAvg = vec3d::Zero();
    for (const auto& subtri : subtris0)
        nhatAvg += 

    for (const auto& subtri : subtris0) {
        // assert(subtri->iVerts[0] == idx0);
        assert(Math::vecEquals(subtri->Xs[0], X0));

        for (int i = 1; i < 3; ++i) {
            const vec3d& rhat = (subtri->Xs[i]-X0).normalized();
            // std::cout << subtri->nhat << ' ' << ehat.cross(rhat) << ' ';
            std::cout << oriented(nhatAvg, ehat, rhat) << '\n';
            // std::cout << acos(ehat.dot(rhat)) << ' ';
        }

        std::cout << '\n';
    }

    std::cout << '\n';
    */
};