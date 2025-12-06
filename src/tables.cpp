#include "tables.h"

void Tables::buildAngularTables(
    const int maxLevel, 
    const double wavenum, 
    const std::vector<realVec>& thetas,
    const std::vector<realVec>& phis) {
    for (int level = 0; level <= maxLevel; ++level) {

        const int nth = thetas[level].size();
        const int nph = 2*nth;

        std::vector<mat3d> ImKK_lvl;
        std::vector<vec3d> kvec_lvl;
        std::vector<mat23d> matToThPh_lvl;
        // std::vector<mat3d> matToThPh_lvl;

        for (int ith = 0; ith < nth; ++ith) {
            const double th = thetas[level][ith];
            for (int iph = 0; iph < nph; ++iph) {
                const double ph = phis[level][iph];

                ImKK_lvl.push_back(Math::IminusRR(th, ph));
                kvec_lvl.push_back(Math::vecSph(wavenum, th, ph));
                matToThPh_lvl.push_back(Math::matToThPh(th, ph));
                // matToThPh_lvl.push_back(Math::matToSph(th, ph));
            }
        }

        ImKK.push_back(ImKK_lvl);
        kvec.push_back(kvec_lvl);
        matToThPh.push_back(matToThPh_lvl);
    }
}

void Tables::buildInterpThetaTable(
    const int maxLevel, const int order, const std::vector<realVec>& thetas) {

    // Do not construct interp table for leaf level
    for (size_t level = 0; level < maxLevel; ++level) { 
        std::vector<realVec> interpTheta_lvl;
        std::vector<size_t> t_lvl;

        const int nth = thetas[level].size();
        const int mth = thetas[level+1].size();

        for (size_t ith = 0; ith < nth; ++ith) {
            realVec interpTheta_lvl_th;
            const double theta = thetas[level][ith];

            // Find idx of child theta nearest parent theta
            const int t = Interp::getNearGLNodeIdx(theta, mth, 0.0, PI);

            // if (level) std::cout << ith << ' ' << theta << ' ' << t << '\n';

            // Assemble child thetas interpolating parent theta
            realVec branchThetas;
            for (int jth = t+1-order; jth <= t+order; ++jth) {

                // Flip jth if not in [0, mth-1]
                size_t jth_flipped = Math::flipIdxToRange(jth, mth);
                auto branchTheta = thetas[level+1][jth_flipped];

                // Extend interpolating thetas to outside [0, \pi] as needed
                if (jth < 0)
                    branchTheta *= -1.0;
                else if (jth >= mth)
                    branchTheta = 2.0*PI - branchTheta;

                branchThetas.push_back(branchTheta);

                // if (level) std::cout << ith << ' ' << jth << ' ' << branchTheta << '\n';
            }
            // if (level) std::cout << '\n';

            for (int k = 0; k <= 2*order-1; ++k) {
                interpTheta_lvl_th.push_back(
                    Interp::evalLagrangeBasis(theta, branchThetas, k));

                // if (level) std::cout << ith << ' ' << k << ' ' << interpTheta_lvl_th[k] << '\n';
            }

            interpTheta_lvl.push_back(interpTheta_lvl_th);
            t_lvl.push_back(t);
        }

        interpTheta.push_back(interpTheta_lvl);

        ts.push_back(t_lvl);
    }
}

void Tables::buildInterpPhiTable(
    const int maxLevel, const int order, const std::vector<realVec>& phis) {

    // Do not construct interp table for leaf level
    for (size_t level = 0; level < maxLevel; ++level) {
        std::vector<realVec> interpPhi_lvl;
        std::vector<size_t> s_lvl;

        const int nph = phis[level].size();
        const int mph = phis[level+1].size();

        for (size_t iph = 0; iph < nph; ++iph) {
            realVec interpPhi_lvl_ph;
            const double phi = phis[level][iph];

            // Find idx of child phi nearest parent phi
            const int s = std::floor(mph * phi / (2.0 * PI));

            // Assemble child phis interpolating parent phi
            realVec branchPhis;
            for (int jph = s+1-order; jph <= s+order; ++jph)
                // Wrap jph if not in [0, mph-1]
                // size_t jph_wrapped = Math::wrapIdxToRange(jph, mph); 
                branchPhis.push_back(2.0*PI*jph/static_cast<double>(mph));

            for (size_t k = 0; k <= 2*order-1; ++k)
                interpPhi_lvl_ph.push_back(
                    Interp::evalTrigBasis(phi, branchPhis, k));
            
            interpPhi_lvl.push_back(interpPhi_lvl_ph);
            s_lvl.push_back(s);
        }

        interpPhi.push_back(interpPhi_lvl);
        ss.push_back(s_lvl);
    }
}

/*
void Tables::buildYlmTables(const int order) {
    auto binom = [](double x, int k) {
        return Math::fallingFactorial(x, k) / Math::factorial(k);
        };

    for (int l = 0; l <= 2*order; ++l) {
        realVec coeffYlm_l, fallingFact_l, legendreSum_l, fracCoeffYlm_l, A_l, Aexp_l;

        for (int m = 0; m <= l; ++m) {
            coeffYlm_l.push_back(Math::coeffYlm(l, m));
            fallingFact_l.push_back(Math::fallingFactorial(l, m));
            legendreSum_l.push_back(binom(l, m) * binom((l+m-1)/2.0, l));
            fracCoeffYlm_l.push_back(sqrt((l-m)/static_cast<double>(l+m)));
        }

        auto pm_l = Math::pm(l);
        for (int m = -l; m <= l; ++m) {
            A_l.push_back(pm_l /
                std::sqrt(static_cast<double>(Math::factorial(l-m)*Math::factorial(l+m))));
            Aexp_l.push_back(1.0 /
                std::sqrt(static_cast<double>(Math::factorial(l-m)*Math::factorial(l+m))));
        }

        coeffYlm_.push_back(coeffYlm_l);
        fallingFact_.push_back(fallingFact_l);
        legendreSum_.push_back(legendreSum_l);
        fracCoeffYlm_.push_back(fracCoeffYlm_l);
        A_.push_back(A_l);
        Aexp_.push_back(Aexp_l);
    }
}*/