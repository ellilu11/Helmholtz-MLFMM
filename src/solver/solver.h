#pragma once

#include "../fmm/fmm.h"
#include "../source.h"

class Solver {

public:
    Solver(const SrcVec&, std::unique_ptr<FMM::Nearfield>);

    virtual void solve(const std::string&) = 0;

    static void printSols(const std::string&, const vecXcd&);

public:
    static vecXcd lvec;
    static vecXcd rvec;
    static vecXcd currents;

protected:
    std::shared_ptr<FMM::Nearfield> nf;
    int numSrcs;
};




