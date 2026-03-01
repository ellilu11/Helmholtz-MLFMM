#include "solver.h"
#include "fmm/nearfield.h"
#include "fmm/node.h"

Solver::Solver(
    SrcVec& srcs,
    std::shared_ptr<FMM::Node> root,
    std::shared_ptr<FMM::Nearfield> nf,
    int maxIter, double EPS)
    : root(std::move(root)),
      nf(std::move(nf)),
      numSrcs(srcs.size()),
      maxIter(maxIter),
      EPS(EPS),
      Qmat(matXcd(numSrcs, 1)),
      gvec(vecXcd::Zero(maxIter+1)),
      vcos(vecXcd::Zero(maxIter)),
      vsin(vecXcd::Zero(maxIter))
{
    /* Sort sources by srcIdx
    std::sort(srcs.begin(), srcs.end(),
        [](std::shared_ptr<Source> src0, std::shared_ptr<Source> src1)
        { return src0->getIdx() < src1->getIdx(); }
    );
    */

    // states.lvec = r = ZI - w = -w assuming I = 0 initially
    // std::transform
    for (int idx = 0; idx < numSrcs; ++idx)
        states.lvec[idx] = -srcs[idx]->getVoltage();

    g0 = states.lvec.norm(); // store g0 for use later
    gvec[0] = g0;

    states.lvec.normalize(); // lvec_0
    Qmat.col(0) = states.lvec; // store lvec as first column of Qmat
}

void Solver::updateRvec(int k) {
    if (!root->isLeaf()) {
        root->mergeMpoleCoeffs();
        root->buildLocalCoeffs();
        FMM::evaluateSols();
    }

    nf->evaluateSols();
}

void Solver::iterateArnoldi(int k) {
    assert(Qmat.cols() == k+1);

    auto& lvec = states.lvec;
    auto& rvec = states.rvec;

    // Do Arnoldi iteration
    vecXcd hcol(k+2);
    for (int i = 0; i <= k; ++i) {
        const auto& q_i = Qmat.col(i); // get lvec_i
        cmplx h = q_i.dot(rvec); // Hermitian dot
        rvec -= h * q_i;
        hcol[i] = h;
    }
    hcol[k+1] = rvec.norm();

    // Replace present lvec with new lvec
    lvec = rvec / hcol[k+1]; // .normalized();

    // Store new lvec as new column of Qmat
    Qmat.conservativeResize(numSrcs, k+2);
    Qmat.col(k+1) = lvec;

    // Store new hcol as new column of Hmat
    Hmat.conservativeResize(k+2, k+1); // 2x1, 3x2, 4x3 ...
    Hmat.row(k+1).setZero();
    Hmat.col(k) = hcol;
}

void Solver::applyGivensRotation(vecXcd& hcol, int k) {
    assert(hcol.size() == k+2);

    for (int i = 0; i <= k-1; ++i) {
        const cmplx hprev = vcos[i] * hcol[i] + vsin[i] * hcol[i+1] ;
        hcol[i+1] = -vsin[i] * hcol[i] + vcos[i] * hcol[i+1];
        hcol[i] = hprev;
    }

    auto [vcos_k, vsin_k] = Math::givensRotation(hcol[k], hcol[k+1]);

    hcol[k] = vcos_k * hcol[k] + vsin_k * hcol[k+1];
    hcol[k+1] = 0.0;

    vcos[k] = vcos_k;
    vsin[k] = vsin_k;
}

void Solver::updateGvec(int k) {

    vecXcd hcol_k = Hmat.col(k);
    applyGivensRotation(hcol_k, k);
    Hmat.col(k) = hcol_k;

    gvec[k+1] = -vsin[k] * gvec[k];
    gvec[k] = vcos[k] * gvec[k];
}

void Solver::printSols(const std::string& fname) {
    namespace fs = std::filesystem;
    fs::path dir = "out/sol";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream file(dir/fname);

    file << std::setprecision(15) << std::scientific;

    // for (const auto& sol : states.rvec) file << sol << '\n';
    for (const auto& curr : states.currents) file << curr << '\n';
}

void Solver::solve(const std::string& fname) {
    auto start = Clock::now();

    int iter = 0;
    do {
        // if (iter && !(iter%10)) std::cout << "   #" << iter << '\n';
        updateRvec(iter);

        iterateArnoldi(iter);

        updateGvec(iter);

        states.rvec = vecXcd::Zero(numSrcs); // reset rvec for next iteration
    } while (abs(gvec[++iter])/g0 > EPS && iter < maxIter); // careful

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n";
    t.printTimes();
    t.resetTimes();

    std::cout << "   # iterations: " << iter << "\n";

    const matXcd& Hp = Hmat.block(0, 0, Hmat.rows()-1, Hmat.cols());
    vecXcd yvec = Hp.lu().solve(gvec.segment(0, iter));
    states.currents = Qmat.leftCols(iter) * yvec;
    std::cout << "   Current norm: " << states.currents.norm() << "\n";

    printSols(fname);
}