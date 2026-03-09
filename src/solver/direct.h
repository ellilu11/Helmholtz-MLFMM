#pragma once

#include "solver.h"

class Direct final : public Solver {

public:
    Direct(const SrcVec& srcs, std::unique_ptr<FMM::Nearfield>);

    void solve(const std::string&) override;
};