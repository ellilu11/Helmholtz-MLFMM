#pragma once

#include "solver.h"

class Direct final : public Solver {

public:
    Direct(SrcVec& srcs,
        std::shared_ptr<FMM::Nearfield>);

    void solve(const std::string&) override;
};