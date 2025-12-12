#include <fstream>
#include "MLFMA.h"
#include "clock.h"
#include "config.h"

using namespace std;

extern auto t = ClockTimes();

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

    // ==================== Build tables ===================== //
    cout << " Building tables...\n";

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

    return 0;

    // ================== Evaluate solutions ================= //
    cout << " Evaluating solutions...\n";
    start = Clock::now();

    Leaf::evaluateSols();

    end = Clock::now();
    duration_ms = end - start;
    Time fmm_duration_ms = end - fmm_start;

    cout << " FMM total elapsed time: " << fmm_duration_ms.count() << " ms\n";

    return 0;
}