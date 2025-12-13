#include <random>
#include "../src/MLFMA.h"
#include "../src/math.h"
#include "../src/node.h"

using namespace std;

shared_ptr<Node> Node::getNode(int nodeIdx) {

    auto node = nodes[nodeIdx];

    cout << "   Selected node at level " << node->getLevel()
         << " with " << node->getIlist().size() << " interaction nodes\n\n";

    return node;
}

void Stem::printLocalCoeffs(std::ofstream& f) {

    for (const auto& coeffs : localCoeffs)
        f << coeffs << '\n';

    /*
    for (const auto& branch : branches)
        branch->printLocalCoeffs(f);
        */
}

void Leaf::printLocalCoeffs(std::ofstream& f) {

    for (const auto& coeffs : localCoeffs)
        f << coeffs << '\n';
}

void testTransl() {
    // ==================== Test translation ===================== //
    cout << " Testing M2L translations...\n";

    realVec testDists =
        { 2, 2.23607, 2.44949, 2.82843, 3, 3.16228, 3.31662, 3.4641,
          3.60555, 3.74166, 4.12311, 4.24264, 4.3589, 4.69042, 5.19615 };

    const auto& translTable = Node::getTables().transl;

    for (const auto& dist : testDists)
        cout << dist << ' ' << translTable[1].at(dist)[0] << '\n';
}

int main() {

    Config config("config/config.txt");

    // ==================== Import geometry ==================== //
    cout << " Constructing RWGs...\n";

    const string nstr = "480";

    auto start = Clock::now();

    auto vertices = importVertices("config/n"+nstr+"/vertices.txt");

    auto tris = importTriangles("config/n"+nstr+"/faces.txt", vertices, config.quadPrec);
    const int numQuads = Triangle::quadPrec2Int(config.quadPrec);

    shared_ptr<Source> Einc = make_shared<Source>(); // initialize incident field

    auto srcs = importRWG("config/n"+nstr+"/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    auto end = Clock::now();
    Time duration_ms = end - start;

    cout << "   # Sources:       " << Nsrcs << '\n';
    cout << "   RWG quad rule:   " << numQuads << "-point\n";
    cout << "   Digit precision: " << config.digits << '\n';
    cout << "   Interp order:    " << config.interpOrder << '\n';
    cout << "   Max node RWGs:   " << config.maxNodeSrcs << '\n';
    cout << fixed << setprecision(3);
    cout << "   Root length:     " << config.rootLeng << '\n';
    cout << "   Wave number:     " << Einc->wavenum << '\n';
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto fmm_start = Clock::now();
    start = Clock::now();

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeSrcs)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    end = Clock::now();
    duration_ms = end - start;

    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Max node level: " << Node::getMaxLvl() << '\n';
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    auto obsNode = root->getNode(55);

    // ==================== Build tables ===================== //
    cout << " Building angular samples...\n";

    start = Clock::now();

    Node::buildAngularSamples();
    Node::buildTables();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Upward pass ===================== //
    cout << " Computing upward pass...\n";

    start = Clock::now();

    root->buildMpoleCoeffs();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Downward pass ==================== //
    cout << " Computing downward pass...\n";
    start = Clock::now();

    root->buildLocalCoeffs();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // Do nearfield test
    ofstream coeffFile("out/lcoeffs.txt");

    obsNode->printLocalCoeffs(coeffFile);

    return 0;

}