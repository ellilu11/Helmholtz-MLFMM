#pragma once

#include "solver.h"

class GMRES final : public Solver {

public:
    GMRES(const SrcVec& srcs,
        std::unique_ptr<FMM::Nearfield>,
        std::shared_ptr<FMM::Node>,
        double, int);

    void updateRvec(int);

    void iterateArnoldi(int);

    void applyGivensRotation(vecXcd&, int);

    void updateGvec(int);

    void solve(const std::string&) override;

private:
    std::shared_ptr<FMM::Node> root;

    // Eigen::SparseLU<sparseMat<cmplx>, Eigen::COLAMDOrdering<int>> precond;
    Eigen::IncompleteLUT<cmplx> precond;

    matXcd Qmat;
    matXcd Hmat;
    vecXcd gvec;
    double g0;

    vecXcd vcos;
    vecXcd vsin;

    double EPS;
    int maxIter;
};



