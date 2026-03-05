#include "direct.h"
#include "../fmm/nearfield.h"
#include "../fmm/node.h"

Direct::Direct(SrcVec& srcs,
    std::shared_ptr<FMM::Nearfield> nf)
    : Solver(srcs, std::move(nf))
{
    // rvec = V = ZI
    for (int idx = 0; idx < numSrcs; ++idx)
        rvec[idx] = srcs[idx]->getVoltage();
};

void Direct::solve(const std::string& fname) {
    auto selfPairs = nf->getSelfPairs();   
    assert(selfPairs.size() == 1); // only one self pair for direct solver

    auto start = Clock::now();

    // Assemble full nearfield matrix
    matXcd Zmat = selfPairs[0].getNearMatrix();
    assert(Zmat.rows() == rvec.rows());
    assert(Zmat.cols() == currents.rows());

    currents = Zmat.lu().solve(rvec);

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n";
    t.printTimes();
    t.resetTimes();

    std::cout << "   Current norm: " << currents.norm() << "\n";

    printSols(fname);
}