#include "../src/interp.h"
#include "../src/leaf.h"
#include "../src/node.h"
#include "../src/stem.h"

using namespace std;

vec3cd Leaf::getLeafSols(const vec3d R) {
    const int nth = thetas[level].size();
    const int nph = phis[level].size();

    vec3cd fld;

    size_t idx = 0;
    for (int ith = 0; ith < nth; ++ith) {
        for (int iph = 0; iph < nph; ++iph) {
            fld += thetaWeights[level][idx] * coeffs[idx];
            idx++;
        }
    }

    return fld;
}

void Node::testFarField(const int nth, const int nph) {

    ofstream obsFile, outFile, outAnlFile
    obsFile.open("config/obss.txt");
    outFile.open("out/ff.txt");
    outAnlFile.open("out/ffAnl.txt");

    const double R = 50.0 * rootLeng;

    for (int ith = 0; ith < nth; ++ith) {
        double th = PI * ith / static_cast<double>(nth);
        for (int iph = 0; iph < nph; ++iph) {
            double ph = 2.0 * PI * iph / static_cast<double>(nph);
            auto obs = vec3d(R, th, ph);
            obss.push_back(obs);
            obsFile << obs << '\n';
        }
    }
}

int main() {
    Config config("config/config.txt");

    // ==================== Import geometry ==================== //
    auto vertices = importVertices("config/n120/vertices.txt");

    auto tris = importTriangles("config/n120/faces.txt", vertices);
    int Ntris;

    auto srcs = importRWG("config/n120/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    cout << " # Sources:           " << Nsrcs << '\n';
    cout << " Root length:         " << config.rootLeng << '\n';
    cout << " Interpolation order: " << config.order << '\n';
    cout << " Max node RWGs:       " << config.maxNodeSrcs << "\n\n";

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
    cout << " Building tables...\n";

    Node::buildAngularSamples();
    Node::buildTables(config);

    // ==================== Test upward pass ===================== //

    root->testFarField(10, 20);

}