#pragma once

#include "srcrwg.h"
#include "subrwg.h"

class BC {
public:
    BC(SrcRWG* const);

private:
    std::vector<SubRWG> subrwgs;
    SrcRWG* const base;

    vec2i numRWGs;
};