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
        [](const std::shared_ptr<Source>& src) { return -src->getVoltage(); });
};

void Solver::printSols(const std::string& fname, const vecXcd& sols) {
    std::filesystem::path dir = "out/sol";
    makeDir(dir);
    std::ofstream file(dir/fname);

    file << std::setprecision(15) << std::scientific;
    for (const auto& sol : sols) file << sol << '\n';
}

void Solver::printScattered(const SrcVec& srcs,
    const std::filesystem::path& dir, const std::string& fname, int nth)
{
    makeDir(dir);
    std::ofstream farfile(dir/fname), thfile(dir/"angles.txt");
    farfile << std::setprecision(15) << std::scientific;

    std::cout << " Computing scattered farfield...\n";

    double k = config.k, k2 = k*k;
    double rcsSum = 0.0;
    for (int ith = 0; ith < nth; ++ith) { // 2*nth to cover great circle
        double theta0 = (ith+0.5)*PI/static_cast<double>(nth); // in [0, 2*pi]
        double theta = (theta0 < PI) ? theta0 : 2*PI - theta0; // fold back to [0, pi]
        double phi = (theta0 < PI) ? 0.0 : PI; // fold back to [0, 2pi]
        assert(theta >= 0.0 && theta <= PI);

        vec3d rhat = Math::fromSph(vec3d(1.0, theta, phi));

        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += currents[src->getIdx()] * src->getFarAlongDir(k*rhat);

        vec2cd far = iu*k*Phys::eta * Math::toThPh(theta, phi) * dirFar;

        double rcs = 4.0*PI/k2 * far.squaredNorm();
        rcsSum += rcs;

        farfile << rcs << '\n'; // squared magnitude of theta and phi components
        thfile << theta0 << ' ' << theta << ' ' << phi << '\n';
    }

    std::cout << "   Mean RCS: "
        << std::setprecision(9) << rcsSum/nth << std::setprecision(3) << "\n";
    std::cout << "   File: " << (dir/fname).generic_string() << "\n\n";
}