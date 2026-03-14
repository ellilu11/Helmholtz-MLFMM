#pragma once

#include "solver.h"

class Direct final : public Solver {

public:
    Direct(const SrcVec& srcs, std::unique_ptr<FMM::Nearfield> nf)
        : Solver(srcs, std::move(nf))
    {
    };

    /* solve()
     * Solve for currents using direct solver, and print solutions to file
     */
    void solve(const std::string& fname) override {
        auto selfPairs = nf->getSelfPairs();
        assert(selfPairs.size() == 1); // only one self pair for direct solver

        std::cout << " Solving for current w/ LU...     ";

        auto start = Clock::now();

        // Convert sparse nearMat to dense matrix for direct solver
        matXcd Zmat = nf->nearMat;
        assert(Zmat.rows() == lvec.rows());
        assert(Zmat.cols() == currents.rows());

        currents = Zmat.lu().solve(lvec);

        Time duration_ms = Clock::now() - start;
        std::cout << " in " << duration_ms.count() << " ms\n";

        std::cout << "   Current norm: "
            << std::setprecision(9) << currents.norm() << std::setprecision(3) << "\n";

        printSols(fname, currents);
    }
};