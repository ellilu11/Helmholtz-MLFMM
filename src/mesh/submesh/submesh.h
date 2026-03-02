#pragma once

namespace Mesh {
    // Constants
    constexpr std::array<double, 5> rcoeffs =
        { -1.0/6.0, 1.0/6.0, -1.0/3.0, 1.0/3.0, 1.0 };

    // Types
    class SrcRWG;
    class SubRWG;
    class BC;

    // Fine mesh data
    std::vector<SubRWG> glSubrwgs;  // list of fine RWGs
    std::vector<double> massCoeffs; // mass coefficients between fine RWGs
    PairHashMap<int> idxMassCoeffs; // subRWG pairs to mass coeff indices

    // Fine mesh maps
    PairHashMap<int> edgeToMid;       // coarse edges to midpoint indices
    PairHashMap<vec2i> fineEdgeToTri; // fine edges to subtri indices 
    PairHashMap<int> fineEdgeToSub;   // fine edges to subrwg indices
    std::vector<std::vector<int>> triToSubs;    // indices of subrwgs containing tri
    std::vector<std::vector<int>> vertToSubs;   // indices of subrwgs containing vert

    // Functions
    void refineMesh(const SrcVec&);
}
