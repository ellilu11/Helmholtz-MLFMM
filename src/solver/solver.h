#pragma once

#include "../fmm/fmm.h"
#include "../source.h"

class Solver {

public:
    Solver(const SrcVec&, std::unique_ptr<FMM::Nearfield>);

    virtual void solve(const std::string&) = 0;

    static void printSols(const std::string&, const vecXcd&);

    static void printScattered(
        const SrcVec& srcs, const std::filesystem::path&, const std::string&, int);

    // static void printSurfCurrents();

public:
    static vecXcd lvec;
    static vecXcd rvec;
    static vecXcd currents;

protected:
    std::unique_ptr<FMM::Nearfield> nf;
    int nsols; // number of unknowns (= nsrcs for PEC)
};




