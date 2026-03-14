#include "solver.h"

vecXcd Solver::currents;
vecXcd Solver::lvec;
vecXcd Solver::rvec;

Solver::Solver(const SrcVec& srcs, std::unique_ptr<FMM::Nearfield> nf)
    : nf(std::move(nf)), numSrcs(srcs.size())
{
    currents = vecXcd::Zero(numSrcs);
    lvec = vecXcd::Zero(numSrcs);
    rvec = vecXcd::Zero(numSrcs);

    // Sort sources by srcIdx to preserve ordering of Zmat, lvec, rvec, and currents
    SrcVec sortedSrcs = srcs;
    std::sort(sortedSrcs.begin(), sortedSrcs.end(),
        [](std::shared_ptr<Source> src0, std::shared_ptr<Source> src1)
        { return src0->getIdx() < src1->getIdx(); }
    );

    // lvec = r = ZI - w = -w assuming I = 0 initially
    std::transform(sortedSrcs.begin(), sortedSrcs.end(), lvec.data(),
        [](const std::shared_ptr<Source>& src) { return -src->getVoltage(); }); // double check sign convention
};

void Solver::printSols(const std::string& fname, const vecXcd& sols) {
    namespace fs = std::filesystem;
    fs::path dir = "out/sol";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream file(dir/fname);

    file << std::setprecision(15) << std::scientific;

    for (const auto& sol : sols) file << sol << '\n';
}