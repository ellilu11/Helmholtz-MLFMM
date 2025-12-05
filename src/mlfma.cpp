#include <fstream>
#include <iostream>
#include "clock.h"
#include "config.h"
#include "MLFMA.h"

using namespace std;

int main() {

    // ===== Tests ===== //
    //for (int l = 1; l <= 5; ++l) {
    //    auto [nodes, weights] = Math::gaussLegendre(l, 1.0E-9, 0.0, PI);
    //    cout << "l = " << l << " ";
    //    for (int k = 0; k < l; ++k)
    //        cout << '(' << nodes[k] << ',' << weights[k] << ") ";
    //    cout << '\n';
    //}
    //return 0;

    // ==================== Import geometry ==================== //
    Config config("config/config.txt");

    auto vertices = importVertices("config/n120/vertices.txt");

    auto tris = importTriangles("config/n120/faces.txt", vertices);
    int Ntris;

    shared_ptr<Src> Einc = make_shared<Src>(); // initialize incident field

    auto srcs = importRWG("config/n120/rwgs.txt", vertices, tris, Einc);
    int Nsrcs = srcs.size();

    Node::setNodeParams(config,Einc);

    cout << " # Sources:           " << Nsrcs << '\n';
    cout << " Root length:         " << config.rootLeng << '\n';
    cout << " Interpolation order: " << config.order << '\n';
    // cout << " Exponential order:   " << Node::getExponentialOrder() << '\n';
    cout << " Max node RWGs:       " << config.maxNodeSrcs << "\n\n";

    // return 0;

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto start = Clock::now();

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeSrcs)
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
    // return 0;

    // ============ Construct directional quantities ========= //
    cout << " Building tables...\n";

    start = Clock::now();

    Node::buildAngularSamples();
    Node::buildTables(config);

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";
    return 0;

    // ==================== Upward pass ===================== //
    cout << " Computing upward pass...\n";

    start = Clock::now();

    root->buildMpoleCoeffs();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    return 0;
}