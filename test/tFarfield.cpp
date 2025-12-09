#include <fstream>
#include <iostream>
#include "../src/MLFMA.h"
#include "../src/interp.h"
#include "../src/leaf.h"
#include "../src/stem.h"

using namespace std;

/* Return observers in Cartesian coordinates */
vector<vec3d> buildObservers(double r, int nth, int nph) {
    ofstream obsFile("config/obss.txt");

    vector<vec3d> obss;
    for (int ith = 0; ith < nth; ++ith) {
        double th = PI * ith / static_cast<double>(nth);
        for (int iph = 0; iph < nph; ++iph) {
            double ph = 2.0 * PI * iph / static_cast<double>(nph);
            auto obs = Math::fromSph(vec3d(r, th, ph));
            obss.push_back(obs);
            obsFile << obs << '\n';
        }
    }

    return obss;
}

vec3cd Leaf::getLeafSol(const vec3d& X) {
    assert(!rwgs.empty());

    const auto [nth, nph] = getNumAngles(level);
    const double phiWeight = 2.0*PI / static_cast<double>(nph);

    const auto dist = X - center;

    // Do angular integral
    vec3cd sol = vec3cd::Zero();
    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        const auto theta = thetas[level][ith];
        const auto thetaWeight = thetaWeights[level][ith];

        vec3cd sol_ph = vec3cd::Zero();
        for (int iph = 0; iph < nph; ++iph) {
            const double phi = phis[level][iph];

            const auto kvec = tables.kvec[level][idx];

            sol_ph += exp(iu*kvec.dot(dist))
                      * tables.matFromThPh[level][idx] 
                      * coeffs[idx];

            idx++;
        }

        sol += thetaWeight * sin(theta) * sol_ph;
    }

    return c0 * mu0 / ( 8.0*PI*PI ) * sol * phiWeight;
}

/*
std::vector<vec2cd> Leaf::getLeafSolsPerTheta(const vec3d& X) {
    assert(!rwgs.empty());

    const auto [nth, nph] = getNumAngles(level);
    const double phiWeight = 2.0*PI / static_cast<double>(nph);

    const vec3d dist = X - center;

    // Do phi integral for each theta
    std::vector<vec2cd> sols;
    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        vec2cd sol_theta = vec2cd::Zero();
        const double theta = thetas[level][ith];

        for (int iph = 0; iph < nph; ++iph) {

            const double phi = phis[level][iph];

            const auto kvec = Math::fromSph(vec3d(wavenum, theta, phi)); // wavenum
            // tables.kvec[level][idx];

            sol_theta +=
                phiWeight
                * exp(iu*kvec.dot(dist))
                * vec2cd(1, 0);// * coeffs[idx]
                ;

            sol_theta[0] +=
                phiWeight * exp(iu*kvec.dot(dist))
                * (coeffs[idx][0] * cos(theta) * cos(phi) - coeffs[idx][1] * sin(phi));

            sol_theta[1] +=
                phiWeight * exp(iu*kvec.dot(dist))
                * (coeffs[idx][0] * cos(theta) * sin(phi) + coeffs[idx][1] * cos(phi));

            idx++;
        }

        sols.push_back(sol_theta);
    }

    return sols;
}*/

/* void Leaf::testFarfieldFromLeaves(const std::vector<vec3d>& obss) {

    const auto [nth, nph] = getNumAngles(0); // get Nth, Nph of root

    string outStr = "out/ff_nl" 
        + to_string(Leaf::getNumLeaves()) 
        + "_nth" + to_string(nth) 
        + "_nph" + to_string(nph)
        + ".txt";

    ofstream outFile(outStr);

    outFile << setprecision(15) << scientific;

    string coeffStr = "out/coeffs_nth" 
        + to_string(nth)
        + "_nph" + to_string(nph)
        + ".txt";

    ofstream coeffFile(coeffStr);

    coeffFile << setprecision(15) << scientific;

    for (const auto& leaf : leaves) {
        leaf->buildMpoleCoeffs();

        for (int idx = 0; idx < nth*nph; ++idx)
            coeffFile << leaf->coeffs[idx] << '\n';

        size_t idx = 0;
        for (int ith = 0; ith < nth; ++ith) {
            for (int iph = 0; iph < nph; ++iph)
                coeffFile << (leaf->coeffs[idx++])[0].real() << ' ';

            coeffFile << '\n';
        }
    }

    for (const auto& obs : obss) {
        vec2cd sol = vec2cd::Zero();

        for (const auto& leaf : leaves) {
            if (leaf->rwgs.empty()) continue;
            sol = sol + leaf->getLeafSol(obs);
        }

        outFile << sol[0].real() << ' ' << sol[1].real() << '\n';

        std::vector<vec2cd> sols;
        for (int i = 0; i < nth; ++i)
            sols.push_back(vec2cd::Zero());

        for (const auto& leaf : leaves) {
            if (leaf->rwgs.empty()) continue;
            sols = sols + leaf->getLeafSolsPerTheta(obs);
        }

        for (int i = 0; i < nth; ++i)
            outFile << (sols[i])[0].real() << ' ';

        outFile << '\n';

    }

    // Also write G-L theta nodes
    ofstream nodeFile("out/nodes_nth" + to_string(nth)+ ".txt");

    for (int lvl = 0; lvl <= 0; ++lvl) {
        for (int ith = 0; ith < nth; ++ith)
            nodeFile << thetas[lvl][ith] << ' ';

        nodeFile << '\n';
    }

}*/

void Leaf::testFarfieldFromLeaves(const std::vector<vec3d>& obss) {

    ofstream outFile("out/ff_nl" + to_string(Leaf::getNumLeaves()) + ".txt");

    outFile << setprecision(15) << scientific;

    for (const auto& leaf : leaves)
        leaf->buildMpoleCoeffs();

    for (const auto& obs : obss) {
        vec3cd sol = vec3cd::Zero();

        for (const auto& leaf : leaves) {
            if (leaf->rwgs.empty()) continue;
            sol = sol + leaf->getLeafSol(obs);
        }

        // Eigen::Array3d solAbs = sol.array().abs();
        vec3d solAbs = sol.cwiseAbs();
        outFile << solAbs << '\n';

        /*outFile << sol[0].real() << ' ' 
                << sol[1].real() << ' '
                << sol[2].real() << '\n';
                */
    }

}

int main() {
    Config config("config/config.txt");

    // ==================== Import geometry ==================== //
    auto vertices = importVertices("config/n120/vertices.txt");

    auto tris = importTriangles("config/n120/faces.txt", vertices);
    int Ntris;

    shared_ptr<Src> Einc = make_shared<Src>(); // initialize incident field

    auto srcs = importRWG("config/n120/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    cout << " # Sources:           " << Nsrcs << '\n';
    cout << " Root length:         " << config.rootLeng << '\n';
    cout << " Interpolation order: " << config.order << '\n';
    cout << " Max node RWGs:       " << config.maxNodeSrcs << "\n";
    cout << " Wave number:         " << Einc->wavenum << "\n\n";

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << Node::getMaxLvl() << '\n';

    // ==================== Build tables ===================== //
    cout << "\n Building tables...\n";

    Node::buildAngularSamples();
    Node::buildTables(config);

    for (int lvl = 0; lvl <= Node::getMaxLvl(); ++lvl) {
        auto [nth, nph] = Node::getNumAngles(lvl);
        cout << " Level " << lvl << ", (Nth,Nph) = "
            << "(" << nth << "," << nph << ")\n";
    }

    // ==================== Test upward pass ===================== //
    cout << "\n Testing upward pass...\n";

    // Set up observers
    const double r = 10.0 * config.rootLeng;
    constexpr int nth = 10;
    constexpr int nph = 20;

    auto obss = buildObservers(r, nth, nph);

    // Do upward pass
    auto start = Clock::now();

    Leaf::testFarfieldFromLeaves(obss);

    // root->testFarField(obss);

    auto end = Clock::now();
    Time duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // Print direct solution
    ofstream outFile("out/ffDir.txt");

    outFile << setprecision(15) << scientific;

    for (const auto& obs : obss) {
        // Eigen::Array3d solAbs = root->getFarSol(obs).array().abs();
        vec3d solAbs = root->getFarSol(obs).cwiseAbs();
        outFile << solAbs << '\n';
    }

    return 0;
}