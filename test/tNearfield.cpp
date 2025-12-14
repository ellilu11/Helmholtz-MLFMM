#include <random>
#include "../src/MLFMA.h"
#include "../src/math.h"
#include "../src/node.h"

using namespace std;

shared_ptr<Node> Node::getNode(int nodeIdx) {

    auto node = nodes[nodeIdx];

    assert(node->isNodeType<Leaf>() && !node->isSrcless());

    cout << "   Selected node at level " << node->getLevel()
         << " with " << node->iList.size() << " interaction nodes and "
         << node->srcs.size() << " RWGs\n\n";

    return node;
}

void Leaf::testMpoleToLocalInLeaf() {

    // Get sols from local coeffs due to iList (assuming L2L is off)
    ofstream outFile("out/nf.txt");

    evalFarSols();

    for (const auto& src : srcs)
        outFile << src->getSol() << '\n';

    // Get sols directly from iList
    ofstream outDirFile("out/nf_dir.txt");

    resetSols();

    for (const auto& node : iList)
        evalPairSols(node);

    for (const auto& src : srcs)
        outDirFile << src->getSol() << '\n';

}

void Node::printAngularSamples(int level) {
    const auto [nth, nph] = getNumAngles(level);

    ofstream thetaFile("out/thetas_nth" + to_string(nth) + ".txt");
    ofstream phiFile("out/phis_nth" + to_string(nth) + ".txt");

    for (int ith = 0; ith < nth; ++ith)
        thetaFile << thetas[level][ith] << '\n';

    for (int iph = 0; iph < nph; ++iph)
        phiFile << phis[level][iph] << '\n';
}


/*void testTransl() {
    // ==================== Test translation ===================== //
    cout << " Testing M2L translations...\n";

    realVec testDists =
        { 2, 2.23607, 2.44949, 2.82843, 3, 3.16228, 3.31662, 3.4641,
          3.60555, 3.74166, 4.12311, 4.24264, 4.3589, 4.69042, 5.19615 };

    const auto& translTable = Node::getTables().transl;

    for (const auto& dist : testDists)
        cout << dist << ' ' << translTable[1].at(dist)[0] << '\n';
}*/

extern auto t = ClockTimes();

int main() {

    Config config("config/config.txt");
    const string configPath = "config/n480/";

    // ==================== Import geometry ==================== //
    auto srcs = importConfig(config, configPath);
    auto Nsrcs = srcs.size();

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

    // ==================== Nearfield test =================== //

    auto obsLevel = obsNode->getLevel();
    auto [nth, nph] = Node::getNumAngles(obsLevel);

    ofstream coeffFile("out/lcoeffs_nth"+to_string(nth)+".txt");
    obsNode->printLocalCoeffs(coeffFile);

    Node::printAngularSamples(obsLevel);

    // auto obsLeaf = dynamic_pointer_cast<Leaf>(obsNode);
    // obsLeaf->testMpoleToLocalInLeaf();

    return 0;

}