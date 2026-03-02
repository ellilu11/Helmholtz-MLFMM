#pragma once

#include "srcrwg.h"
#include "submesh.h"
#include "subrwg.h"

class Mesh::BC {
public:
    BC(SrcRWG* const);

    void receiveRvals();

private:
    SrcRWG* const base;           // pointer to base srcRWG
    std::array<std::vector<int>,2> iSubss;  // indices of subRWGs at both common vertices
    std::array<std::vector<double>, 2> pcoeffss; // precomputed fine RWG -> BC coeffs

    cmplx rval;
};