#include "gmres.h"
#include "../fmm/nearfield.h"
#include "../fmm/node.h"

GMRES::GMRES(
    const SrcVec& srcs,
    std::unique_ptr<FMM::Nearfield> nf,
    std::shared_ptr<FMM::Node> root,
    double EPS, int maxIter)
    : Solver(srcs, std::move(nf)),
      root(std::move(root)),
      Qmat(matXcd(numSrcs, 1)),
      gvec(vecXcd::Zero(maxIter+1)),
      vcos(vecXcd::Zero(maxIter)),
      vsin(vecXcd::Zero(maxIter)),
      EPS(EPS),
      maxIter(maxIter)
{
    currents = vecXcd::Zero(numSrcs);
    lvec = vecXcd::Zero(numSrcs);
    rvec = vecXcd::Zero(numSrcs);

    // Sort sources by srcIdx
    SrcVec sortedSrcs = srcs;
    std::sort(sortedSrcs.begin(), sortedSrcs.end(),
        [](std::shared_ptr<Source> src0, std::shared_ptr<Source> src1)
        { return src0->getIdx() < src1->getIdx(); }
    );

    auto M = this->nf->getNearMatrix(numSrcs);
    // precond.compute(M); // compute LU of M for preconditioning
    precond.analyzePattern(M);
    precond.factorize(M);

    // lvec = r = ZI - w = -w assuming I = 0 initially
    // std::transform
    for (int idx = 0; idx < numSrcs; ++idx)
        lvec[idx] = -sortedSrcs[idx]->getVoltage();
    lvec = precond.solve(lvec); // apply M^{-1} here if preconditioning is desired

    g0 = lvec.norm(); // store g0 for use later
    gvec[0] = g0;
 
    lvec.normalize(); // lvec_0

    Qmat.col(0) = lvec; // store lvec as first column of Qmat
}

void GMRES::updateRvec(int k) {
    if (!root->isLeaf()) {
        root->mergeMpoleCoeffs();
        root->buildLocalCoeffs();
        FMM::evaluateSols();
    }
    nf->evaluateSols();

    // Apply M^{-1} here if preconditioning is desired
    rvec = precond.solve(rvec);
}

void GMRES::iterateArnoldi(int k) {
    assert(Qmat.cols() == k+1);

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

void GMRES::applyGivensRotation(vecXcd& hcol, int k) {
    assert(hcol.size() == k+2);

    for (int i = 0; i <= k-1; ++i) {
        cmplx hprev = vcos[i] * hcol[i] + vsin[i] * hcol[i+1] ;
        hcol[i+1] = -vsin[i] * hcol[i] + vcos[i] * hcol[i+1];
        hcol[i] = hprev;
    }

    auto [vcos_k, vsin_k] = Math::givensRotation(hcol[k], hcol[k+1]);

    hcol[k] = vcos_k * hcol[k] + vsin_k * hcol[k+1];
    hcol[k+1] = 0.0;

    vcos[k] = vcos_k;
    vsin[k] = vsin_k;
}

void GMRES::updateGvec(int k) {
    vecXcd hcol_k = Hmat.col(k);
    applyGivensRotation(hcol_k, k);
    Hmat.col(k) = hcol_k;

    gvec[k+1] = -vsin[k] * gvec[k];
    gvec[k] = vcos[k] * gvec[k];
}

void GMRES::solve(const std::string& fname) {
    // If maxIter = 0, just do one MVM and exit
    if (!maxIter) {
        updateRvec(0);
        printSols(fname, rvec);
        std::cout << '\n';
        return;
    }

    std::string method = root->isLeaf() ? "direct... " : "FMM...    ";
    std::cout << " Solving for current w/ " << method;

    auto start = Clock::now();
    int iter = 0;
    do {
        if (iter && !(iter%100)) std::cout << " #" << iter << ' ';
        updateRvec(iter);
        iterateArnoldi(iter);
        updateGvec(iter);
        rvec = vecXcd::Zero(numSrcs); // reset rvec for next iteration
    } while (abs(gvec[++iter])/g0 > EPS && iter < maxIter); // careful

    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n";
    t.printTimes();
    t.resetTimes();

    std::cout << "   # Iterations: " << iter << "\n";

    const matXcd& Hp = Hmat.block(0, 0, Hmat.rows()-1, Hmat.cols());
    vecXcd yvec = Hp.lu().solve(gvec.segment(0, iter));
    currents = Qmat.leftCols(iter) * yvec;
    std::cout << "   Current norm: " 
        << std::setprecision(9) << currents.norm() << std::setprecision(3) << "\n";

    printSols(fname, currents);
}