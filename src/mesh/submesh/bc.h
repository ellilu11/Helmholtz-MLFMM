#pragma once

#include "../srcrwg.h"
#include "subrwg.h"

class Mesh::BC {
public:
    BC(SrcRWG* const);

    void receiveRvals();

private:
    SrcRWG* const base;           // pointer to base srcRWG
    std::array<intVec,2> iSubss;  // indices of subRWGs at both common vertices
    std::array<realVec, 2> pcoeffss; // precomputed fine RWG -> BC coeffs

    cmplx rval;
};