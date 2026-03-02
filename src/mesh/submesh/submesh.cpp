#include "submesh.h"

void Mesh::refineMesh(const SrcVec& rwgs) {
    Triangle::refineVertices();
    Triangle::refineTriangles();
    Triangle::buildEdgeToTri();

    SrcRWG::refineRWGs();

    for (const auto& rwg : rwgs)
        dynamic_pointer_cast<SrcRWG>(rwg)->findSubRWGs();

    SubRWG::buildVertsToSubRWGs(nverts);

    for (const auto& rwg : rwgs)
        dynamic_pointer_cast<SrcRWG>(rwg)->buildBC();

    SubRWG::buildMassCoeffs();

    edgeToMid.clear();
    fineEdgeToTri.clear();
    fineEdgeToSub.clear();
    vertToSubs.clear();
}