#include "direct.h"
#include "../fmm/nearfield.h"
#include "../fmm/node.h"

Direct::Direct(const SrcVec& srcs,
    std::shared_ptr<FMM::Nearfield> nf)
    : Solver(srcs, std::move(nf))
{
    currents = vecXcd::Zero(numSrcs);
    rvec = vecXcd::Zero(numSrcs);

    // rvec = V = ZI
    for (int idx = 0; idx < numSrcs; ++idx)
        rvec[idx] = -srcs[idx]->getVoltage(); // Double check sign convention here
};

void Direct::solve(const std::string& fname) {
    auto selfPairs = nf->getSelfPairs();
    assert(selfPairs.size() == 1); // only one self pair for direct solver

    std::cout << " Solving for current w/ direct LU...";

    auto start = Clock::now();

    matXcd Zmat = selfPairs[0].getNearMatrix();
    assert(Zmat.rows() == rvec.rows());
    assert(Zmat.cols() == currents.rows());

    currents = Zmat.lu().solve(rvec);

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n";

    std::cout << "   Current norm: "
        << std::setprecision(9) << currents.norm() << std::setprecision(3) << "\n";

    printSols(fname, currents);
}