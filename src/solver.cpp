#include "solver.h"

Solver::Solver(const SrcVec& srcs, std::shared_ptr<Node> root, int numIter)
    : root(std::move(root)),
      numSols(srcs.size()),
      numIter(numIter),
      Q(matXcd(numSols, 1)),
      gvec(vecXcd::Zero(numIter+1)),
      currents(vecXcd::Zero(numSols)), // Guess zero currents
      qvec(std::make_shared<vecXcd>(vecXcd::Zero(numSols))),
      sols(std::make_shared<vecXcd>(vecXcd::Zero(numSols)))
{
    // Guess zero currents
    // std::transform
    for (int idx = 0; idx < numSols; ++idx)
        (*qvec)[idx] = -srcs[idx]->getVoltage();

    gvec[0] = (*qvec).norm();

    (*qvec).normalize(); // qvec_0
    Q.col(0) = *qvec; // store qvec as first column of Q
}

void Solver::updateSols(int k) {

    std::cout << " Do FMM iteration #" << k << '\n';

    auto iter_start = Clock::now();

    root->buildMpoleCoeffs();

    root->buildLocalCoeffs();

    Leaf::evaluateSols();

    Time fmm_duration_ms = Clock::now() - iter_start;
    if (!k) {
        std::cout << "   Elapsed time (S2M): " << t.S2M.count() << " ms\n";
        std::cout << "   Elapsed time (M2M): " << t.M2M.count() << " ms\n";
        std::cout << "   Elapsed time (M2L): " << t.M2L.count() << " ms\n";
        std::cout << "   Elapsed time (L2L): " << t.L2L.count() << " ms\n";
    }
    std::cout << "   Total elapsed time: " << fmm_duration_ms.count() << " ms\n\n";
}

void Solver::iterateArnoldi(int k) {

    assert(Q.cols() == k+1);

    vecXcd hvec(k+2);

    // Do Arnoldi iteration
    for (int i = 0; i <= k; ++i) {
        const auto& q_i = Q.col(i); // get qvec_i
        cmplx h = q_i.dot(*sols); // Hermitian dot
        *sols -= h * q_i;
        hvec[i] = h;
    }
    hvec[k+1] = (*sols).norm();

    // Replace present qvec with new qvec
    *qvec = (*sols).normalized();

    // Store new qvec as new column of Q
    Q.conservativeResize(numSols, k+2);
    Q.col(k+1) = *qvec;

    // Store new hvec as new column of H
    H.conservativeResize(k+2, k+1); // 2x1, 3x2, 4x3 ...
    H.col(k) = hvec;
}

void Solver::applyGivensRotation(
    vecXcd& hvec, vecXcd& vcos, vecXcd& vsin, int k) {

    assert(hvec.size() == k+2);

    for (int i = 0; i < k; ++i) {
        const cmplx htemp = vcos[i] * hvec[i] + vsin[i] * hvec[i+1] ;
        hvec[i+1] = -vsin[i] * hvec[i] + vcos[i] * hvec[i+1];
        hvec[i] = htemp;
    }

    auto [vcos_k, vsin_k] = Math::givensRotation(hvec[k], hvec[k+1]);

    hvec[k] = vcos_k * hvec[k] + vsin_k * hvec[k+1];
    hvec[k+1] = 0.0;

    vcos[k] = vcos_k;
    vsin[k] = vsin_k;
}

void Solver::updateGvec(vecXcd& vcos, vecXcd& vsin, int k) {

    vecXcd H_k = H.col(k);
    applyGivensRotation(H_k, vcos, vsin, k);
    H.col(k) = H_k;

    gvec[k+1] = -vsin[k] * gvec[k];
    gvec[k] = vsin[k] * gvec[k];
}

void Solver::solve() {

    //
    namespace fs = std::filesystem;
    fs::path dir = "out/sol";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream file(dir/"sol.txt", std::ios::app);
    //

    vecXcd vcos = vecXcd::Zero(numIter);
    vecXcd vsin = vecXcd::Zero(numIter);

    for (int iter = 0; iter < numIter; ++iter) {
        updateSols(iter);

        iterateArnoldi(iter);

        updateGvec(vcos, vsin, iter);

        resetSols();

        // 
        const auto& Hp = H.block(0, 0, H.rows()-1, H.cols());
        vecXcd yvec = Hp.lu().solve(gvec.segment(0, iter+1));

        currents = Q.leftCols(iter+1) * yvec;

        printSols(file);
        //
    }
}

// void Solver::printSols(const std::string& fname)
void Solver::printSols(std::ofstream& os)
{

    // os << std::setprecision(9) << std::scientific;

    //for (const auto& sol : *sols)
    //    file << sol.real() << ' ' << sol.imag() << '\n';

    for (const auto& curr : currents)
        // file << curr.real() << ' ' << curr.imag() << '\n';
        os << curr.real() << ' ';
    os << '\n';
}