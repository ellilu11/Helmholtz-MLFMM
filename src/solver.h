#pragma once

#include <filesystem>
#include "fmm/fmm.h"
#include "source.h"

class Solver {

public :
    Solver(
        SrcVec& srcs, 
        std::shared_ptr<FMM::Node>, 
        int, double);

    void updateRvec(int);

    void iterateArnoldi(int);

    void applyGivensRotation(vecXcd&, int);

    void updateGvec(int);

    void solve(const std::string&);

    void printSols(const std::string&);

private :
    std::shared_ptr<FMM::Node> root;

    int numSrcs;
    int maxIter;
    double EPS;

    matXcd Qmat;
    matXcd Hmat;
    vecXcd gvec;
    double g0;

    vecXcd vcos;
    vecXcd vsin;
};