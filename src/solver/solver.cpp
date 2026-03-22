#include "solver.h"

vecXcd Solver::currents;
vecXcd Solver::lvec;
vecXcd Solver::rvec;

Solver::Solver(const SrcVec& srcs, std::unique_ptr<FMM::Nearfield> nf)
    : nf(std::move(nf)), nsols(srcs.size())
{
    currents = vecXcd::Zero(nsols);
    lvec = vecXcd::Zero(nsols);
    rvec = vecXcd::Zero(nsols);

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

/* printScattered()
 * Compute and print the scattered far field in the direction of a great circle
 * centered at the origin, with angles theta in [0, pi] and phi = 0 or pi
 * (fold back to [0, pi] to get full coverage of great circle)
 */
void Solver::printScattered(const SrcVec& srcs,
    const std::filesystem::path& dir, const std::string& fname, int nangles)
{
    makeDir(dir);
    std::ofstream farfile(dir/fname), thfile(dir/"angles.txt");
    farfile << std::setprecision(15) << std::scientific;

    std::cout << " Computing scattered farfield...\n";

    double k = config.k, k2 = k*k;
    double rcsSum = 0.0;
    for (int ith = 0; ith < 2*nangles; ++ith) { // 2*nth to cover great circle
        double theta0 = (ith+0.5)*PI/static_cast<double>(nangles); // in [0, 2*pi]
        double theta = (theta0 < PI) ? theta0 : 2*PI - theta0; // fold back to [0, pi]
        double phi = (theta0 < PI) ? 0.0 : PI; // fold back to [0, 2pi]
        assert(theta >= 0.0 && theta <= PI);
    //for (int iph = 0; iph < nangles; ++iph) {
    //    double theta = PI / 2.0;
    //    double phi = (iph+0.5)*PI/static_cast<double>(nangles); // in [0, pi]

        vec3d rhat = Math::fromSph(vec3d(1.0, theta, phi));

        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += currents[src->getIdx()] * src->getFarAlongDir(k*rhat);

        // Get theta and phi components of scattered far field
        vec2cd far = iu*k*Phys::eta * Math::toThPh(theta, phi) * dirFar;
        double rcs = 4.0*PI/k2 * far.squaredNorm();

        // Get V or H component of scattered far field
        //cmplx far = iu*k*Phys::eta * dirFar[0];
        //double rcs = 4.0*PI/k2 * std::norm(far); 

        rcsSum += rcs;

        farfile << rcs << '\n'; // squared magnitude of theta and phi components
        thfile << theta0 << ' ' << theta << ' ' << phi << '\n';
        // thfile << phi << ' ' << theta << ' ' << phi << '\n';
    }

    std::cout << "   Mean RCS: "
        << std::setprecision(9) << rcsSum/nangles << std::setprecision(3) << "\n";
    std::cout << "   File: " << (dir/fname).generic_string() << "\n\n";
}

/* printSurfCurrents()
 * Get the surface current on this triangle by summing contributions from adjacent RWGs
 * Use the center of the triangle as the evaluation point for the RWG function
cmplx Solver::printSurfCurrents() const {
    auto triToRWG = triToRWGs[iTri];

    cmplx J = 0.0;
    for (int i = 0; i < 3; ++i) {
        int iRWG = triToRWG.iRWGs[i];
        int isMinus = triToRWG.isMinus[i];
        auto rwg = glSrcs[iRWG];
        const auto& Xnc = rwg->getVertsNC();

        double rwgFunc =
            Math::sign(isMinus) * rwg->leng / (2.0*area) * (center - Xnc[isMinus]);

        J += Solver::currents[iRWG];
    }
    return J;
}
*/