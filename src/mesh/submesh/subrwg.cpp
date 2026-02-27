#include "subrwg.h"

Mesh::SubRWG::SubRWG(int iSub, const vec4i& idx4)
    : RWG(nullptr, iSub, idx4)
{
    const auto& [tri0, tri1] = getTris();

    // If an RWG straddles two coarse mesh vertices, 
    // it does not contribute to the BC function at either vertex
    if (tri0.iVerts[0] == tri1.iVerts[0])
        iVertsCoarse = tri0.iVerts[0];

    // Record triangles contained in this subRWG
    triToSubs[tri0.iTri].push_back(iSub);
    triToSubs[tri1.iTri].push_back(iSub);

    /*std::cout << "Built subRWG #" << iSrc << " w/ common vertices # "
        << iVertsC[0] << ' '<< iVertsC[1] << " and non-common vertices # "
        << iVertsNC[0] << ' ' << iVertsNC[1] << " and coarse vertex # "
        << (iVertsCoarse.has_value() ? std::to_string(iVertsCoarse.value()) : "N/A") << "\n";*/

    /* Reorder tris if needed
    const vec3d& nhat0 = dX.cross(Xnc[0] - Xc[0]);
    const vec3d& nhat1 = dX.cross(Xnc[1] - Xc[0]);
    assert(nhat0.dot(nhat0 - nhat1) > 0);
    // std::cout << nhat0.dot(nhat0 - nhat1) << '\n';
    // if (nhat0.dot(nhat0 - nhat1) < 0) std::swap(tris[0],tris[1]);
    */
}

void Mesh::SubRWG::buildVertsToSubRWGs(int numVerts) {
    // For each subrwg, map its coarse mesh vertex to itself
    vertToSubs.resize(numVerts);

    int iSub = 0;
    // Find subRWGs contributing to BC at coarse mesh vertices
    for (const auto& rwg : glSubrwgs) {

        if (rwg.iVertsCoarse)
            vertToSubs[rwg.iVertsCoarse.value()].push_back(iSub);
        ++iSub;
    }

    /*
    int vIdx = 0;
    for (const auto& vertToSub : vertToSubs) {
        std::cout << "Vertex #" << vIdx++ << " has subRWGs #";
        for (const auto& iSub : vertToSub) {
            std::cout << ' ' << iSub;
        }
        std::cout << " with common vertices # ";
        for (const auto& iSub : vertToSub) {
            const auto& rwg = glSubrwgs[iSub];
            std::cout << '(' << rwg.iVertsC[0] << ',' << rwg.iVertsC[1] << ") ";
        }

        std::cout << '\n';
    }
    */
}

void Mesh::SubRWG::buildMassCoeffs() {
    massCoeffs.reserve(glSubrwgs.size()*glSubrwgs.size());

    size_t idx = 0;
    for (const auto& rwg : glSubrwgs) {
        // Compute self coeffs
        massCoeffs.push_back(rwg.getMassCoeff(rwg));
        idxMassCoeffs.emplace(
            makeUnordered(rwg.iSrc, rwg.iSrc), idx++);

        // Compute pairwise coeffs
        for (auto iTri : rwg.iTris) {
            // std::cout << "Triangle #" << iTri << " is contained inside RWGs: ";

            // Only loop over RWGs sharing a triangle with this RWG
            for (auto iSubs : triToSubs[iTri]) {
                const auto& other = glSubrwgs[iSubs];

                if (rwg.iSrc < other.iSrc) {
                    massCoeffs.push_back(rwg.getMassCoeff(other));
                    idxMassCoeffs.emplace(
                        makeUnordered(rwg.iSrc, other.iSrc), idx++);
                }
            }
        }
    }
    assert(massCoeffs.size() == idx);
    
    massCoeffs.shrink_to_fit();
    triToSubs.clear();

    std::cout << std::setprecision(6);
    for (int i = 0; i < glSubrwgs.size(); ++i) {
        for (int j = 0; j < glSubrwgs.size(); ++j) {
            const auto key = makeUnordered(i,j);
            if (idxMassCoeffs.find(key) != idxMassCoeffs.end()) {
                std::cout << massCoeffs[idxMassCoeffs[key]] << ' ';
            } else {
                std::cout << "0.000000 ";
            }
        }
        std::cout << '\n';
    }
}

double Mesh::SubRWG::getMassCoeff(const SubRWG& other) const {
    const auto& Xnc0 = getVertsNC();
    const auto& Xnc1 = other.getVertsNC();
    
    // double massNum = 0.0;
    double mass = 0.0;

    size_t iPair0 = 0;
    for (const auto& iTri0 : iTris) {
        const auto& tri0 = glTris[iTri0];

        size_t iPair1 = 0;
        for (const auto& iTri1 : other.iTris) {
            if (iTri0 == iTri1) {
                // Numeric
                //for (const auto& [node, weight] : tri0.getQuads())
                //    massNum += weight * (node - Xnc0[iPair0]).dot(node - Xnc1[iPair1])
                //            * Math::sign(iPair0) * Math::sign(iPair1)
                //            / (2.0 * tri0.area);

                // Analytic
                const auto& [X0, X1, X2] = tri0.getVerts();
                const double sum =
                    (X0.squaredNorm() + X1.squaredNorm() + X2.squaredNorm() + 
                     X0.dot(X1) + X0.dot(X2) + X1.dot(X2)) / 12.0
                    - (Xnc0[iPair0] + Xnc1[iPair1]).dot(tri0.center) / 2.0
                    + Xnc0[iPair0].dot(Xnc1[iPair1]) / 2.0;

                mass += sum * Math::sign(iPair0) * Math::sign(iPair1)
                            / (2.0 * tri0.area);

            }

            ++iPair1;
        }

        ++iPair0;
    }

    // assert(Math::approxEquals(massNum, mass));
    // std::cout << "Mass coeff between rwg #" << iSrc << " and rwg #" << other.iSrc << ": " << mass << '\n';

    return leng * other.leng * mass;
}

void Mesh::SubRWG::setOriented(int iVert, const vec3d& nhat, const vec3d& ehat) {
    const vec3d& Xref = glVerts[iVert]; // vertex on coarse mesh
    const auto& [X0,X1] = getVertsC(); // vertices of common edge

    auto [i0, i1] = iVertsC;
    assert(i0 == iVert || i1 == iVert);

    const vec3d& rhat = (i0 == iVert ? X1-X0 : X0-X1).normalized();
    const double angle = atan2(nhat.dot(ehat.cross(rhat)),ehat.dot(rhat));
    oriented = (Math::fless(angle,0.0) ? angle+2.0*PI : angle);
}