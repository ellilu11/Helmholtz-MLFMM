#include "level.h"

std::vector<double> FMM::Level::dists;
std::vector<vec3d> FMM::Level::rhats;
std::array<vec3d, 316> FMM::Level::dXs;

void FMM::Level::buildAngularSamples()
{
    double nodeLeng = Mesh::rootLeng / pow(2.0, level);

    // Use excess bandwidth formula
    int tau = ceil((1.73*config.k*nodeLeng +
        2.16*pow(config.digits, 2.0/3.0)*pow(config.k*nodeLeng, 1.0/3.0)));

    L = floor(0.50*tau); // TODO: Find optimal formula for M2L series truncation number

    // Construct thetas
    int nth = tau+1;
    std::tie(thetas, weights) = Math::gaussLegendre(nth, 0.0, PI);

    // Absorb sin(theta) into weights
    std::transform(weights.begin(), weights.end(), thetas.begin(), weights.begin(),
        [](double weight, double theta) { return weight * sin(theta); }
    );

    // Construct phis
    int nph = 2*nth;
    phis.resize(nph);

    for (int iph = 0; iph < nph; ++iph)
        phis[iph] = 2.0*PI*iph/static_cast<double>(nph);

    // std::cout << "   (" << level << "," << thetas.size() << "," << phis.size() << ")\n";
}

void FMM::Level::buildAngularMatrices() {
    auto [nth, nph] = getNumAngles();
    size_t nDirs = nth*nph;

    khat.resize(nDirs);
    toThPh.resize(nDirs);
    ImRR.resize(nDirs);

    size_t iDir = 0;
    for (int ith = 0; ith < nth; ++ith) {
        double theta = thetas[ith];

        for (int iph = 0; iph < nph; ++iph) {
            double phi = phis[iph];

            khat[iDir] = Math::fromSph(vec3d(1.0, theta, phi));
            toThPh[iDir] = Math::toThPh(theta, phi);
            ImRR[iDir] = Math::ImRR(khat[iDir]);

            ++iDir;
        }
    }
}

std::vector<interpPair> FMM::Level::getInterpTheta(int srcLvl, int tgtLvl)
{
    int order = config.interpOrder;

    const std::vector<double>&
        srcThetas = fmmLevels[srcLvl].thetas, tgtThetas = fmmLevels[tgtLvl].thetas;
    int mth = srcThetas.size(), nth = tgtThetas.size();

    std::vector<interpPair> interpPairs;
    interpPairs.reserve(nth);

    for (size_t jth = 0; jth < nth; ++jth) {
        double tgtTheta = tgtThetas[jth];

        int nearIdx = Math::getNearGLNodeIdx(tgtTheta, mth, 0.0, PI);

        // Assemble source thetas interpolating target theta
        std::vector<double> interpThetas(2*order);
        for (int ith = nearIdx+1-order, k = 0; ith <= nearIdx+order; ++ith, ++k) {

            // Flip ith if not in [0, mth-1]
            int ith_flipped = Math::flipIdxToRange(ith, mth);

            double srcTheta = srcThetas[ith_flipped];

            // Extend source thetas to outside [0, pi] as needed
            if (ith < 0) srcTheta *= -1.0;
            else if (ith >= mth) srcTheta = 2.0*PI - srcTheta;

            interpThetas[k] = srcTheta;
        }

        vecXd coeffs(2*order);
        for (int k = 0; k < 2*order; ++k)
            coeffs[k] = Math::evalLagrangeBasis(tgtTheta, interpThetas, k);

        interpPairs.emplace_back(coeffs, nearIdx);
    }

    return interpPairs;
}

std::vector<interpPair> FMM::Level::getInterpPhi(int srcLvl, int tgtLvl)
{
    int order = config.interpOrder;

    const std::vector<double>&
        srcPhis = fmmLevels[srcLvl].phis, tgtPhis = fmmLevels[tgtLvl].phis;
    int mph = srcPhis.size(), nph = tgtPhis.size();

    std::vector<interpPair> interpPairs;
    interpPairs.reserve(nph);

    for (size_t jph = 0; jph < nph; ++jph) {
        double tgtPhi = tgtPhis[jph];

        int nearIdx = std::floor(mph * tgtPhi / (2.0*PI));

        // Assemble source phis interpolating target phi
        std::vector<double> interpPhis(2*order);
        for (int iph = nearIdx+1-order, k = 0; iph <= nearIdx+order; ++iph, ++k)
            interpPhis[k] = 2.0*PI*iph/static_cast<double>(mph);

        vecXd coeffs(2*order);
        for (int k = 0; k < 2*order; ++k)
            coeffs[k] = Math::evalLagrangeBasis(tgtPhi, interpPhis, k);

        interpPairs.emplace_back(coeffs, nearIdx);
    }

    return interpPairs;
}

Map<vecXcd> FMM::Level::getAlpha() {
    using namespace Math;

    int nth = thetas.size();
    int nps = std::floor(config.overInterp*(nth-1));
    double nodeLeng = config.rootLeng / pow(2.0, level);

    Map<vecXcd> alpha;
    for (const auto& dist : dists) {
        double kr = k * dist * nodeLeng;

        vecXcd transl_dist(nps);
        for (int ips = 0; ips < nps; ++ips) {
            double xi = cos(PI*ips/static_cast<double>(nps-1));
            cmplx coeff = 0.0;

            for (int l = 0; l <= L; ++l)
                coeff += powI(l) * (2.0*l+1.0)
                * sphericalHankel1(kr, l)
                * legendreP(xi, l).first;

            transl_dist[ips] = iu * k / (4.0*PI) * coeff;
        }

        alpha.emplace(dist, transl_dist);
    }

    assert(alpha.size() == dists.size());
    return alpha;
};

HashMap<interpPair> FMM::Level::getInterpPsi() {
    int order = config.interpOrder;

    // Find all unique psi = acos(khat.dot(rhat))
    const auto [nth, nph] = getNumAngles();
    int nDir = nth*nph;

    std::vector<double> psis(nDir*rhats.size());

    size_t m = 0;
    for (size_t iDir = 0; iDir < nDir; ++iDir) {
        const auto& khat = this->khat[iDir];

        for (const auto& rhat : rhats)
            psis[m++] = acos(khat.dot(rhat));
    }

    std::sort(psis.begin(), psis.end());
    psis.erase(std::unique(psis.begin(), psis.end()), psis.end());

    // Compute Lagrange coefficients for each possible psi
    int nps = std::floor(config.overInterp*(nth-1));

    HashMap<interpPair> interpPairs;
    for (auto psi : psis) {
        // Find idx of psi node nearest this psi
        int nearIdx = std::floor((nps-1) * psi / PI);

        // Assemble psis interpolating this psi
        std::vector<double> psis(2*order);
        for (int ips = nearIdx+1-order, k = 0; ips <= nearIdx+order; ++ips, ++k)
            psis[k] = PI*ips/static_cast<double>(nps-1);

        // CONSIDER: Barycentric coordinates
        vecXd coeffs(2*order);
        for (size_t k = 0; k < 2*order; ++k)
            coeffs[k] = Math::evalLagrangeBasis(psi, psis, k);

        interpPairs.emplace(psi, std::make_pair(coeffs, nearIdx));
    }

    // assert(interpPairs.size() == psis.size());

    return interpPairs;
}

void FMM::Level::buildTranslationTable() {
    int order = config.interpOrder;

    const auto& alphas = getAlpha();
    const auto& interpPsis = getInterpPsi();

    const auto [nth, nph] = getNumAngles();
    int nDir = nth*nph;

    int nps = std::floor(config.overInterp*(nth-1));

    transl.reserve(dXs.size());
    for (const auto& dX : dXs) {
        double r = dX.norm();
        const auto& rhat = dX / r;

        const auto& alpha_dX = alphas.at(r);

        arrXcd transl_dX(nDir);
        for (int iDir = 0; iDir < nDir; ++iDir) {
            const auto& khat = this->khat[iDir];
            double psi = acos(khat.dot(rhat));
            const auto [interpPsi, nearIdx] = interpPsis.at(psi);

            cmplx translCoeff = 0.0;
            for (int ips = nearIdx+1-order, k = 0; k < 2*order; ++ips, ++k) {
                int ips_flipped = Math::flipIdxToRange(ips, nps);

                translCoeff += alpha_dX[ips_flipped] * interpPsi[k];
            }

            transl_dX(iDir) = translCoeff;
        }

        transl.emplace(dX, transl_dX);
    }

    assert(transl.size() == dXs.size());
}
