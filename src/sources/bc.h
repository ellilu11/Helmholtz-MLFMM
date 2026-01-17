#pragma once

#include "srcrwg.h"
#include "subrwg.h"

class BC {
public:
    BC(SrcRWG* const);

private:
    // std::array<std::shared_ptr<RWG>, 10> subrwgs;
    SubRWGVec subrwgs;
    SrcRWG* const base;

    vec2i numRWGs;
};