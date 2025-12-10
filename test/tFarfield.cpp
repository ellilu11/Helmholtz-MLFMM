#include <fstream>
#include <iostream>
#include "../src/MLFMA.h"
#include "../src/interp.h"
#include "../src/leaf.h"
#include "../src/stem.h"

using namespace std;

/* Build observers on sphere at sampled directions */
/*vector<vec3d> Node::getObssAtAngularSamples(double r) {

    const auto [nth, nph] = getNumAngles(level);

    vector<vec3d> obss;
    for (int ith = 0; ith < nth; ++ith) {
        const auto theta = thetas[level][ith];

        for (int iph = 0; iph < nph; ++iph) {
            const double phi = phis[level][iph];

            auto obs = Math::fromSph(vec3d(r, theta, phi));

            obss.push_back(obs);
        }
    }
}*/

std::vector<vec3cd> Leaf::getLeafSols(double r) {
    assert(!rwgs.empty());

    const cmplx C = -iu * c0 * wavenum * mu0
        * exp(iu*wavenum*r) / (4.0*PI*r);

    const auto [nth, nph] = getNumAngles(level);

    std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const auto kvec = tables.kvec[level][idx];
            
            sols[idx] = C * exp(-iu*kvec.dot(center)) * coeffs[idx];

            idx++;
        }
    }

    return sols;
}

/* testFarfieldFromLeaves(r)
 * Print total farfield along leaf sampled directions at distance r,
 * assuming all non-empty leaves are at same level
 */ 
void Leaf::testFarfieldFromLeaves(double r) {

    ofstream outFile("out/ff_lvl" + to_string(maxLevel) + ".txt");

    outFile << setprecision(15) << scientific;

    const auto [nth, nph] = getNumAngles(maxLevel);

    std::vector<vec3cd> sols(nth*nph, vec3cd::Zero());

    for (const auto& leaf : leaves) {
        if (leaf->rwgs.empty()) continue;
        
        // only works if all non-empty leaves are at maxlevel
        assert(leaf->level == maxLevel); 

        leaf->buildMpoleCoeffs();

        sols = sols + leaf->getLeafSols(r);
    }

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const vec3d solAbs = sols[idx].cwiseAbs();

            outFile << solAbs << '\n';

            idx++;
        }
    }
}

void Node::testFarfieldDir(double r) {

    ofstream outFile("out/ffDir.txt");

    outFile << setprecision(15) << scientific;

    const auto [nth, nph] = getNumAngles(level);

    auto sols = getFarSols(r);

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {

        for (int iph = 0; iph < nph; ++iph) {

            const vec3d solAbs = sols[idx].cwiseAbs();

            outFile << solAbs << '\n';

            idx++;
        }
    }
}

void Node::printAngularSamples() {
    ofstream thetaFile("out/thetas_lvl" + to_string(maxLevel) + ".txt");
    ofstream phiFile("out/phis_lvl" + to_string(maxLevel) + ".txt");

    const auto [nth, nph] = getNumAngles(maxLevel);

    for (int ith = 0; ith < nth; ++ith)
        thetaFile << thetas[maxLevel][ith] << '\n';

    for (int iph = 0; iph < nph; ++iph)
        phiFile << phis[maxLevel][iph] << '\n';
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

    const int maxLevel = Node::getMaxLvl();

    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        auto [nth, nph] = Node::getNumAngles(lvl);
        cout << "   (Lvl,Nth,Nph) = "
            << "(" << lvl << "," << nth << "," << nph << ")\n";
    }

    // ==================== Test upward pass ===================== //
    cout << "\n Testing upward pass...\n";

    const double r = 20.0*config.rootLeng; // pick r >> rootLeng

    // Compute farfield from multipole coeffs
    auto start = Clock::now();

    Leaf::testFarfieldFromLeaves(r);
    // root->testFarField(r);

    auto end = Clock::now();
    Time duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // Compute farfield directly
    // root->testFarfieldDir(r);

    // Print out theta and phi samples at max level
    Node::printAngularSamples();

    return 0;
}