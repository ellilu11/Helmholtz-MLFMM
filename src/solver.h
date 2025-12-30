#pragma once

#include <filesystem>
#include "types.h"
#include "sources/source.h"

class Node;

class Solver {

public :
    Solver(const SrcVec& srcs, std::shared_ptr<Node>, int);

    void iterateArnoldi(int);

    void updateSols(int);

    void applyGivensRotation(vecXcd&, vecXcd&, vecXcd&, int);

    void updateGvec(vecXcd&, vecXcd&, int);

    void solve();

    std::shared_ptr<vecXcd> getQvec() { return qvec; }

    std::shared_ptr<vecXcd> getSols() { return sols; }

    // cmplx getQvec(size_t idx) { return qvec[idx]; }

    // void addToSols(size_t idx, cmplx val) { sols[idx] += val; }

    void resetSols() { (*sols) = vecXcd::Zero(numSols); }

    // void printSols(const std::string&);
    void printSols(std::ofstream&);

private :
    //Solver();
    //Solver(const Solver&) = delete;
    //Solver& operator=(const Solver&) = delete;
    
    std::shared_ptr<Node> root;

    int numSols;
    int numIter;

    matXcd Q;
    matXcd H;

    vecXcd gvec;
    vecXcd currents;

    std::shared_ptr<vecXcd> qvec;
    std::shared_ptr<vecXcd> sols;
    
};