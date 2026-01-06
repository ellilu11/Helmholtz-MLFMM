#include "../src/MLFMA.h"
#include "../src/fileio.h"
#include "../src/math.h"
#include "../src/node.h"

using namespace std;

extern auto t = ClockTimes();

cmplx sphFunc(double th, double ph) {
    // return 0.3989422804 * exp(iu*ph);

    // return 0.48860251190292 * cos(th);

    return 0.345494149471336 * sin(th) * exp(iu*ph);
};

/*void Stem::tInterpPhi(int level) {
    const int order = config.interpOrder;

    // Evaluate function at source nodes
    const auto mph = getNumAngles(level+1).second;
    cmplxVec vals;

    for (int iph = 0; iph < mph; ++iph) {
        const double phi = phis[level+1][iph];
        vals.push_back(phiFunc(phi));
    }

    // Interpolate over phi
    const auto nph = getNumAngles(level).second;
    cmplxVec interpVals(nph, 0.0);
    size_t n = 0;

    for (int jph = 0; jph < nph; ++jph) {
        const auto [interp, nearIdx] = tables.interpPhi[level][jph];

        for (int iph = nearIdx+1-order, k = 0; iph <= nearIdx+order; ++iph, ++k) {

            const int iph_wrapped = Math::wrapIdxToRange(iph, mph);

            interpVals[n] += interp[k] * vals[iph_wrapped];
        }

        ++n;
    }

    // Do inner product (weighted)
    const double phiWeight = 2.0*PI / static_cast<double>(nph);
    cmplx intVal = 0.0;

    for (int jph = 0; jph < nph; ++jph) {
        const double phi = phis[level][jph];
        intVal += phiWeight * conj(phiFunc(phi)) * interpVals[jph];
    }

    std::cout << "Integrated val using interp: " << std::setprecision(15) << intVal << '\n';
}

void Stem::tAnterpPhi(int level) {
    const int order = config.interpOrder;

    // Evaluate function times weights at source nodes
    const auto mph = getNumAngles(level).second;
    cmplxVec vals;

    const double phiWeight = 2.0*PI / static_cast<double>(mph);

    for (int iph = 0; iph < mph; ++iph) {
        const double phi = phis[level][iph];

        vals.push_back(phiWeight*phiFunc(phi));
    }

    // Anterpolate over phi
    const auto nph = getNumAngles(level+1).second;
    cmplxVec anterpVals(nph, 0.0);
    cmplxVec anterpLongVals(nph+2*order, 0.0);

    for (int iph = 0; iph < mph; ++iph) { // over parent phis anterpolating child phis

        const auto [interp, nearIdx] = tables.interpPhi[level][iph];

        // if (!ith) std::cout << nearIdx << ' ';

        for (int jph = -order; jph < nph+order; ++jph) { // over child phis to anterpolate

            // shift from jph \in [nearIdx+1-order,nearIdx+order] to k \in [0,2*order-1]
            const int k = jph - (nearIdx+1-order);

            // if iph \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            // if (!ith) std::cout << k << ' ';

            anterpLongVals[jph+order] += interp[k] * vals[iph];
        }

        for (int jph = 0; jph < nph; ++jph) {

            anterpVals[jph] = anterpLongVals[jph+order];

            if (jph < order)
                anterpVals[jph] += anterpLongVals[jph+order+nph];
            else if (jph >= nph-order)
                anterpVals[jph] += anterpLongVals[jph+order-nph];
        }
    }

    // Do inner product
    cmplx intVal = 0.0;

    for (int jph = 0; jph < nph; ++jph) {
        const double phi = phis[level+1][jph];
        intVal += conj(anterpVals[jph]) * phiFunc(phi);
    }

    std::cout << "Integrated val using anterp: " << std::setprecision(15) << intVal << '\n';
}*/

void Stem::tInterp(int level) {
    const int order = config.interpOrder;

    // Evaluate function at source nodes
    const auto [mth, mph] = getNumAngles(level+1);
    cmplxVec vals;

    for (int ith = 0; ith < mth; ++ith) {
        const double theta = thetas[level+1][ith];

        for (int iph = 0; iph < mph; ++iph) {
            const double phi = phis[level+1][iph];

            vals.push_back(sphFunc(theta, phi));
        }
    }

    // Interpolate function values to target nodes
    const auto [nth, nph] = getNumAngles(level);
    cmplxVec innerVals(nth*mph, 0.0);

    size_t m = 0;
    for (int jth = 0; jth < nth; ++jth) {
        const auto [interp, nearIdx] = tables.interpTheta[level][jth];

        for (int iph = 0; iph < mph; ++iph) {

            for (int ith = nearIdx+1-order, k = 0; ith <= nearIdx+order; ++ith, ++k) {

                // Flip ith if not in [0, mth-1]
                const int ith_flipped = Math::flipIdxToRange(ith, mth);

                const bool outOfRange = ith != ith_flipped; // jth < 0 || jth >= mth;

                int iph_shifted = iph;

                // if theta \notin [0, pi] then if:
                // phi \in (0, pi) add pi, phi \in (pi, 2pi) subtract pi
                if (outOfRange)
                    iph_shifted += ((iph < mph/2) ? mph/2 : -mph/2);

                const int m_shifted = ith_flipped*mph + iph_shifted;

                // if (outOfRange)
                // if (ith < ith_flipped)
                innerVals[m] += interp[k] * vals[m_shifted];
            }

            ++m;
        }
    }

    // Interpolate over phi
    cmplxVec interpedVals(nth*nph, 0.0);
    size_t n = 0;
    for (int jth = 0; jth < nth; ++jth) {
        const int jthmph = jth*mph;

        for (int jph = 0; jph < nph; ++jph) {
            const auto [interp, nearIdx] = tables.interpPhi[level][jph]; // TODO: don't need to lookup for every jth

            for (int iph = nearIdx+1-order, k = 0; iph <= nearIdx+order; ++iph, ++k) {

                // Wrap iph if not in [0, mph-1]
                const int iph_wrapped = Math::wrapIdxToRange(iph, mph);

                interpedVals[n] += interp[k] * innerVals[jthmph + iph_wrapped];
            }

            ++n;
        }
    }

    // Do inner product (weighted)
    const double phiWeight = 2.0*PI / static_cast<double>(nph);
    cmplx intVal = 0.0;
    size_t l = 0;
    for (int jth = 0; jth < nth; ++jth) {
        const double theta = thetas[level][jth];
        const double thetaWeight = thetaWeights[level][jth];

        for (int jph = 0; jph < nph; ++jph) {
            const double phi = phis[level][jph];

            intVal += thetaWeight * phiWeight * conj(sphFunc(theta, phi)) * interpedVals[l++];
        }
    }
    assert(l == interpedVals.size());

    std::cout << "Integrated val using interp: " << std::setprecision(15) << intVal << '\n';
}

void Stem::tAnterp(int level) {
    const int order = config.interpOrder;

    // Evaluate function times weights at source nodes
    const auto [mth, mph] = getNumAngles(level);
    cmplxVec vals;

    const double phiWeight = 2.0*PI / static_cast<double>(mph);
    for (int ith = 0; ith < mth; ++ith) {
        const double theta = thetas[level][ith];
        const double thetaWeight = thetaWeights[level][ith];

        for (int iph = 0; iph < mph; ++iph) {
            const double phi = phis[level][iph];

            vals.push_back(phiWeight*thetaWeight*sphFunc(theta, phi));
        }
    }

    // Anterpolate function values to target nodes on extended grid
    const auto [nth, nph] = getNumAngles(level+1);
    assert(!(nph%2));

    const int Nth = nth+2*order, Nph = nph+2*order;

    // Anterpolate over phi
    cmplxVec innerVals(mth*Nph, 0.0);

    for (int iph = 0; iph < mph; ++iph) { // over parent phis anterpolating child phis
        const auto [interp, nearIdx] = tables.interpPhi[level][iph];

        for (int jph = -order; jph < nph+order; ++jph) { // over child phis to anterpolate

            // shift from jph \in [nearIdx+1-order,nearIdx+order] to k \in [0,2*order-1]
            const int k = jph - (nearIdx+1-order);

            // if iph \notin [nearIdx+1-order,nearIdx+order], matrix element is zero
            if (k < 0 || k >= 2*order) continue;

            for (int ith = 0; ith < mth; ++ith) // over parent thetas (unanterpolated)
                innerVals[ith*Nph+jph+order] += interp[k] * vals[ith*mph+iph];
        }
    }

    // Anterpolate over theta
    cmplxVec longVals(Nth*Nph, 0.0);

    for (int ith = 0; ith < mth; ++ith) { // over parent thetas anterpolating child thetas

        const auto [interp, nearIdx] = tables.interpTheta[level][ith];

        for (int jph = -order; jph < nph+order; ++jph) { // over child phis (anterpolated)

            for (int jth = -order; jth < nth+order; ++jth) { // over child thetas to anterpolate

                // shift from ith \in [t+1-order,t+order] to k \in [0,2*order-1]   
                const int k = jth - (nearIdx+1-order);

                // if ith \notin [t+1-order,t+order], matrix element is zero
                if (k < 0 || k >= 2*order) continue;

                longVals[(jth+order)*Nph+jph+order] += interp[k] * innerVals[ith*Nph+jph+order];
            }
        }
    }

    // Contract grid points in theta and phi
    cmplxVec anterpVals(nth*nph, 0.0);
    int dirIdx = 0;

    for (int jth = 0; jth < nth; ++jth) {
        for (int jph = 0; jph < nph; ++jph) {
            anterpVals[dirIdx] = longVals[(jth+order)*Nph+jph+order];

            // Handle nodes near prime meridian
            int jph_wrapped = jph;
            if (jph < order) {
                jph_wrapped += nph;
                anterpVals[dirIdx] += longVals[(jth+order)*Nph+jph_wrapped+order];
            } else if (jph >= nph-order) {
                jph_wrapped -= nph;
                anterpVals[dirIdx] += longVals[(jth+order)*Nph+jph_wrapped+order];
            }

            // Handle nodes near poles
            int jph_shifted = jph;
            jph_shifted += ((jph < nph/2) ? nph/2 : -nph/2);
            assert(jph_shifted >=0 && jph_shifted < nph);

            if (jth < order) {
                int jth_flipped = -jth-1;
                assert(jth_flipped+order >= 0);
                anterpVals[dirIdx] += longVals[(jth_flipped+order)*Nph+jph_shifted+order];
            }
            else if (jth >= nth-order) {
                int jth_flipped = 2*nth-jth-1;
                anterpVals[dirIdx] += longVals[(jth_flipped+order)*Nph+jph_shifted+order];
            }


            ++dirIdx;
        }
    }


    /*
    for (int jph = 0; jph < nph; ++jph) {
        const int ithnph = ith*nph;

        innerVals[ithnph+jph] = innerLongVals[ithnph+jph+order];

        if (jph < order)
            innerVals[ithnph+jph] += innerLongVals[ithnph+jph+order+nph];
        else if (jph >= nph-order)
            innerVals[ithnph+jph] += innerLongVals[ithnph+jph+order-nph];
    }*/

    // Do inner product
    cmplx intVal = 0.0;
    size_t l = 0;
    for (int jth = 0; jth < nth; ++jth) {
        const double theta = thetas[level+1][jth];

        for (int jph = 0; jph < nph; ++jph) {
            const double phi = phis[level+1][jph];

            intVal += conj(anterpVals[l++]) * sphFunc(theta, phi);
        }
    }
    assert(l == anterpVals.size());

    std::cout << "Integrated val using anterp: " << std::setprecision(15) << intVal << '\n';
}

int main() {
    // ===================== Read config ==================== //
    Config config("config/config.txt");

    auto [srcs, Einc] = importFromConfig(config);
    auto nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto fmm_start = Clock::now();
    auto start = Clock::now();

    shared_ptr<Node> root;
    if (nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    auto end = Clock::now();
    Time duration_ms = end - start;

    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << Node::getMaxLvl() << '\n';
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build tables ===================== //
    cout << " Building angular samples...\n";

    start = Clock::now();

    Node::buildAngularSamples();
    Node::buildTables();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // Do test
    const int lvl = 0;
    Stem::tInterp(lvl);
    Stem::tAnterp(lvl);

    return 0;
}

/*void Stem::testInvInterp(int level) {

    // Evaluate function at source nodes
    const auto [mth, mph] = getNumAngles(level);
    cmplxVec vals;

    for (int ith = 0; ith < mth; ++ith) {
        const double theta = thetas[level][ith];

        for (int iph = 0; iph < mph; ++iph) {
            const double phi = phis[level][iph];

            vals.push_back(sphFunc(theta, phi));
        }
    }

    // Interpolate function values to target nodes
    const auto [nth, nph] = getNumAngles(level+1);
    const int order = config.interpOrder;
    cmplxVec innerVals(nth*mph, 0.0);

    size_t m = 0;
    for (int jth = 0; jth < nth; ++jth) {
        const int t = tables.invIdxTheta[level][jth];

        for (int iph = 0; iph < mph; ++iph) {

            for (int ith = t+1-order, k = 0; ith <= t+order; ++ith, ++k) {

                const int ith_flipped = Math::flipIdxToRange(ith, mth);

                const bool outOfRange = ith != ith_flipped;

                int iph_shifted = iph;

                if (outOfRange)
                    iph_shifted += ((iph < mph/2) ? mph/2 : -mph/2);

                const int m_shifted = ith_flipped*mph + iph_shifted;

                innerVals[m] +=
                    tables.invInterpTheta[level][jth][k] * vals[m_shifted];
            }

            m++;
        }
    }

    cmplxVec interpedVals(nth*nph, 0.0);
    size_t n = 0;
    for (int jth = 0; jth < nth; ++jth) {

        for (int jph = 0; jph < nph; ++jph) {
            const int s = tables.invIdxPhi[level][jph];

            for (int iph = s+1-order, k = 0; iph <= s+order; ++iph, ++k) {

                const int iph_wrapped = Math::wrapIdxToRange(iph, mph);

                interpedVals[n] +=
                    tables.invInterpPhi[level][jph][k]
                    * innerVals[jth*mph + iph_wrapped];
            }

            n++;
        }
    }

    ofstream outFile0("out/valsInterped.txt");
    for (int idx = 0; idx < nth*nph; ++idx)
        outFile0 << interpedVals[idx].real() << ' ' << interpedVals[idx].imag() << '\n';

    // Evaluate function at target nodes directly
    ofstream outFile1("out/valsDirect.txt");

    for (int jth = 0; jth < nth; ++jth) {
        const double theta = thetas[level+1][jth];

        for (int jph = 0; jph < nph; ++jph) {
            const double phi = phis[level+1][jph];

            const cmplx val = sphFunc(theta, phi);

            outFile1 << val.real() << ' ' << val.imag() << '\n';
        }
    }

}*/