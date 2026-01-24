#pragma once

#include "srcrwg.h"
#include "subrwg.h"

class Mesh::BC {
public:
    BC(SrcRWG* const);

    void receiveRvals();

private:
    SrcRWG* const base;
    std::array<intVec,2> iSubss;
    cmplx rval;
};