#pragma once

#include "solver.h"

class GMRES final : public Solver {

public:
    GMRES(SrcVec& srcs,
        std::shared_ptr<FMM::Nearfield>,
        std::shared_ptr<FMM::Node>,
        double, int);

    void updateRvec(int);

    void iterateArnoldi(int);

    void applyGivensRotation(vecXcd&, int);

    void updateGvec(int);

    void solve(const std::string&) override;

private:
    std::shared_ptr<FMM::Node> root;

    matXcd Qmat;
    matXcd Hmat;
    vecXcd gvec;
    double g0;

    vecXcd vcos;
    vecXcd vsin;

    double EPS;
    int maxIter;
};



