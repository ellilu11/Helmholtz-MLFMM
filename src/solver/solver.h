#pragma once

#include "../fmm/fmm.h"
#include "../source.h"

class Solver {

public:
    Solver(SrcVec&,
        std::shared_ptr<FMM::Nearfield>);

    virtual void solve(const std::string&) = 0;

    void printSols(const std::string&);

public:
    static vecXcd lvec;
    static vecXcd rvec;
    static vecXcd currents;

protected:
    std::shared_ptr<FMM::Nearfield> nf;

    int numSrcs;
};




