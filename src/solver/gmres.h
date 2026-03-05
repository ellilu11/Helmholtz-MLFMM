#pragma once

#include "solver.h"

class GMRES final : public Solver {

public:
    GMRES(SrcVec& srcs,
        std::shared_ptr<FMM::Nearfield>,
        std::shared_ptr<FMM::Node>,
        int, double);

    void updateRvec(int);

    void iterateArnoldi(int);

    void applyGivensRotation(vecXcd&, int);

    void updateGvec(int);

    void solve(const std::string&) override;

private:
    std::shared_ptr<FMM::Node> root;

    int maxIter;
    double EPS;

    matXcd Qmat;
    matXcd Hmat;
    vecXcd gvec;
    double g0;

    vecXcd vcos;
    vecXcd vsin;
};



