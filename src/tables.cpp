#include "tables.h"

void Tables::buildAngularTables() {
   
    for (int level = 0; level <= Node::maxLevel; ++level) {

        const auto [nth, nph] = Node::getNumAngles(level);

        std::vector<mat3d> ImKK_lvl;
        std::vector<vec3d> khat_lvl;
        std::vector<mat23d> matToThPh_lvl;
        std::vector<mat32d> matFromThPh_lvl;
        
        // std::vector<mat3d> matToSph_lvl;
        // std::vector<mat3d> matFromSph_lvl;

        for (int ith = 0; ith < nth; ++ith) {
            const double th = Node::thetas[level][ith];

            for (int iph = 0; iph < nph; ++iph) {
                const double ph = Node::phis[level][iph];

                const auto& khat = Math::fromSph(vec3d(1.0,th,ph));

                khat_lvl.push_back(khat);
                ImKK_lvl.push_back(Math::IminusRR(khat));
                
                matToThPh_lvl.push_back(Math::matToThPh(th, ph));
                matFromThPh_lvl.push_back(Math::matFromThPh(th, ph));

                // matToSph_lvl.push_back(Math::matToSph(th, ph));
                // matFromSph_lvl.push_back(Math::matFromSph(th, ph));

                // if (level == 2 && ith == 0) 
                    // std::cout << std::setprecision(9) << th << ' ' << ph << ' ' << khat << "\n\n";
                //    std::cout << std::setprecision(9) << Math::IminusRR(khat) << "\n\n";
            }
        }

        ImKK.push_back(ImKK_lvl);
        khat.push_back(khat_lvl);
        matToThPh.push_back(matToThPh_lvl);
        matFromThPh.push_back(matFromThPh_lvl);

        // matToSph.push_back(matToSph_lvl);
        // matFromSph.push_back(matFromSph_lvl);
    }
}

void Tables::buildInterpThetaTable() {

    for (size_t level = 0; level < Node::maxLevel; ++level) { 
        std::vector<realVec> interpTheta_lvl;
        std::vector<int> t_lvl;

        const int nth = Node::getNumAngles(level).first;
        const int mth = Node::getNumAngles(level+1).first;

        for (size_t jth = 0; jth < nth; ++jth) {
            realVec interpTheta_lvl_th;
            const double theta = Node::thetas[level][jth];

            const int t = Interp::getNearGLNodeIdx(theta, mth, 0.0, PI);

            // Assemble child thetas interpolating parent theta
            // TODO: Use splicing with std::span
            realVec branchThetas;
            for (int ith = t+1-order; ith <= t+order; ++ith) {

                // Flip ith if not in [0, mth-1]
                int ith_flipped = Math::flipIdxToRange(ith, mth);

                auto branchTheta = Node::thetas[level+1][ith_flipped];

                // Extend interpolating thetas to outside [0, pi] as needed
                if (ith < 0)
                    branchTheta *= -1.0;
                else if (ith >= mth)
                    branchTheta = 2.0*PI - branchTheta;

                branchThetas.push_back(branchTheta);
            }

            for (int k = 0; k < 2*order; ++k)
                interpTheta_lvl_th.push_back(
                    Interp::evalLagrangeBasis(theta, branchThetas, k));

            interpTheta_lvl.push_back(interpTheta_lvl_th);
            t_lvl.push_back(t);
        }

        interpTheta.push_back(interpTheta_lvl);

        ts.push_back(t_lvl);
    }
}

void Tables::buildInterpPhiTable() {

    for (size_t level = 0; level < Node::maxLevel; ++level) {
        std::vector<realVec> interpPhi_lvl;
        std::vector<int> s_lvl;

        const int nph = Node::getNumAngles(level).second;
        const int mph = Node::getNumAngles(level+1).second;

        for (size_t jph = 0; jph < nph; ++jph) {
            realVec interpPhi_lvl_ph;
            const double phi = Node::phis[level][jph];

            const int s = std::floor(mph * phi / (2.0 * PI));

            // Assemble child phis interpolating parent phi
            // TODO: Use splicing with std::span
            realVec branchPhis;
            for (int iph = s+1-order; iph <= s+order; ++iph) {

                auto branchPhi = 2.0*PI*iph/static_cast<double>(mph);

                branchPhis.push_back(branchPhi);
            }

            for (size_t k = 0; k < 2*order; ++k)
                interpPhi_lvl_ph.push_back(
                    Interp::evalLagrangeBasis(phi, branchPhis, k));
            
            interpPhi_lvl.push_back(interpPhi_lvl_ph);
            s_lvl.push_back(s);
        }

        interpPhi.push_back(interpPhi_lvl);
        ss.push_back(s_lvl);
    }
}

void Tables::buildTranslationTable() {

    using namespace Math;

    const int rootLeng = Node::config.rootLeng;
    const double wavenum = Node::wavenum;

    const auto& dists = getINodeDistances();

    for (size_t level = 0; level <= Node::maxLevel; ++level) {

        const int L = Node::Ls[level];

        const int nth = Node::getNumAngles(level).first;
        const int nps = std::floor(Q*(nth-1));

        const double nodeLeng = rootLeng / pow(2.0, level);

        Map<double,vecXcd> transl_lvl;
        
        for (const auto& dist : dists) {

            const double kr = wavenum * dist * nodeLeng;

            vecXcd transl_dist(nps);

            for (int ips = 0; ips < nps; ++ips) {

                const double psi = PI*ips/static_cast<double>(nps-1);
                const double xi = cos(psi);
                // const double xi = 2.0*ips/static_cast<double>(nps-1)-1.0;
                // const double xi = -2.0*ips/static_cast<double>(nps-1)+1.0;
                cmplx coeff = 0.0;

                for (int l = 0; l <= L; ++l)
                    coeff += 
                        powI(l) * (2.0*l+1.0)
                        * sphericalHankel1(kr, l) 
                        * legendreP(xi, l).first;

                transl_dist[ips] = iu * wavenum / (4.0*PI) * coeff;
            }

            transl_lvl.emplace(dist, transl_dist);
        }

        transl.push_back(transl_lvl);
    }

};


void Tables::buildInterpPsiTable() { // CONSIDER: Interp over xi = cos(psi)

    //std::cout << "   Finding all psi\n";
    //auto start = Clock::now();

    const auto& rhats = Math::getINodeDirections();

    std::vector<realVec> psis;

    // Find all unique psi = acos(khat.dot(rhat)) at each level
    for (size_t level = 0; level <= Node::maxLevel; ++level) {

        const auto [nth, nph] = Node::getNumAngles(level);

        realVec psis_lvl;

        size_t idx = 0;
        for (size_t ith = 0; ith < nth; ++ith) {
            for (size_t iph = 0; iph < nph; ++iph) {

                const auto& khat_ = khat[level][idx++];

                // Loop over all possible rhat
                for (const auto& rhat : rhats)
                    psis_lvl.push_back(acos(khat_.dot(rhat)));
                    // psis_lvl.push_back(khat_.dot(rhat));

            }
        }

        std::sort(psis_lvl.begin(), psis_lvl.end()); 

        psis_lvl.erase(
            std::unique(psis_lvl.begin(), psis_lvl.end()), psis_lvl.end());

        // for (auto psi : psis_lvl) std::cout << psi << ' ';

        psis.push_back(psis_lvl);

    }

    //auto end = Clock::now();
    //Time duration_ms = end - start;
    //std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    //std::cout << "   Computing Lagrange coeffs for each psi\n";
    //start = Clock::now();

    // Compute Lagrange coefficients for each possible psi at each level
    for (size_t level = 0; level <= Node::maxLevel; ++level) {

        const int nth = Node::getNumAngles(level).first;
        const int nps = std::floor(Q*(nth-1));

        HashMap<double,vecXcd> interpPsi_lvl;
        HashMap<double,int> s_lvl;

        size_t idx = 0;
        for (auto psi : psis[level]) {

            // Find idx of psi node nearest this psi
            const int s = std::floor((nps-1) * psi / PI);
            // const int s = std::floor((nps-1) * (psi + 1.0) / 2.0);
            // const int s = std::floor(-(nps-1) * (psi - 1.0) / 2.0);

            // Assemble psis interpolating this psi
            // TODO: Use splicing with std::span
            realVec psis;
            for (int ips = s+1-order; ips <= s+order; ++ips)
                psis.push_back(PI*ips/static_cast<double>(nps-1));
                // psis.push_back(2.0*ips/static_cast<double>(nps-1)-1.0);
                // psis.push_back(-2.0*ips/static_cast<double>(nps-1)+1.0);

            vecXd interpPsi_ps(2*order);

            for (size_t k = 0; k < 2*order; ++k)
                interpPsi_ps[k] = 
                     Interp::evalLagrangeBasis(psi, psis, k);

            interpPsi_lvl.emplace(psi, interpPsi_ps);
            s_lvl.emplace(psi, s);

        }

        interpPsi.push_back(interpPsi_lvl);
        ssps.push_back(s_lvl);
    }

    //end = Clock::now();
    //duration_ms = end - start;
    //std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

}
