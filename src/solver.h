#pragma once

#include <filesystem>
#include "fmm/fmm.h"
#include "source.h"

class Solver {

public :
    Solver(
        SrcVec& srcs, 
        std::shared_ptr<FMM::Node>,
        std::shared_ptr<FMM::Nearfield>,
        int, double);

    void updateRvec(int);

    void iterateArnoldi(int);

    void applyGivensRotation(vecXcd&, int);

    void updateGvec(int);

    void solve(const std::string&);

    void printSols(const std::string&);

private :
    std::shared_ptr<FMM::Node> root;
    std::shared_ptr<FMM::Nearfield> nf;

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