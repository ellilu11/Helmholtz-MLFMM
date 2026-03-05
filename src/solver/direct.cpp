#include "direct.h"
#include "../fmm/nearfield.h"
#include "../fmm/node.h"

Direct::Direct(SrcVec& srcs,
    std::shared_ptr<FMM::Nearfield> nf)
    : Solver(srcs, std::move(nf))
{
    currents = vecXcd::Zero(numSrcs);
    rvec = vecXcd::Zero(numSrcs);

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

    // Print out Zmat
    std::ofstream zfile("out/test/zmat.txt");
    for (int i = 0; i < Zmat.rows(); ++i) {
        for (int j = 0; j < Zmat.cols(); ++j)
            zfile << Zmat(i, j) << ' ';
        zfile << '\n';
    }

    // Print out rvec
    std::ofstream rfile("out/test/rvec.txt");
    for (int i = 0; i < rvec.size(); ++i)
        rfile << rvec[i] << '\n';

    currents = Zmat.lu().solve(rvec);

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n";

    std::cout << "   Current norm: " << currents.norm() << "\n";

    // std::cout << std::setprecision(9) << (Zmat * currents - rvec).transpose() << "\n";

    printSols(fname);
}