#include <random>
#include "../src/MLFMA.h"
#include "../src/import.h"
#include "../src/math.h"
#include "../src/node.h"

using namespace std;

shared_ptr<Node> Node::getNode() {

    size_t nodeIdx = 0;

    auto node = nodes[nodeIdx];

    while (node->isNodeType<Stem>() || node->iList.empty() || node->isSrcless())
        node = nodes[++nodeIdx];

    // assert(node->isNodeType<Leaf>() && !node->isSrcless());

    cout << "   Selected node at level " << node->getLevel()
         << " with " << node->iList.size() << " interaction nodes and "
         << node->srcs.size() << " srcs\n\n";

    return node;
}

void Leaf::testMpoleToLocalInLeaf() {

    auto [nth, nph] = Node::getNumAngles(level);

    // Get sols from local coeffs due to iList (assuming L2L is off)
    ofstream outFile("out/nf/nf_nth"+to_string(nth)+".txt");

    evalFarSols();

    for (const auto& src : srcs)
        outFile << src->getSol() << '\n';

    // Get sols directly from iList
    ofstream outDirFile("out/nf/nf_dir.txt");

    resetSols();

    for (const auto& node : iList)
        evalPairSols(node);

    for (const auto& src : srcs)
        outDirFile << src->getSol() << '\n';

}

void Node::printAngularSamples(int level) {
    const auto [nth, nph] = getNumAngles(level);

    ofstream thetaFile("out/nf/thetas_nth" + to_string(nth) + ".txt");
    ofstream phiFile("out/nf/phis_nth" + to_string(nth) + ".txt");

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
    // ===================== Read config ==================== //
    Config config("config/config.txt");

    auto [srcs, Einc] = importFromConfig(config);
    auto Nsrcs = srcs.size();

    Node::setNodeParams(config, Einc);

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto fmm_start = Clock::now();
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

    auto obsNode = root->getNode();

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

    //ofstream coeffFile("out/nf/lcoeffs_nth"+to_string(nth)+".txt");
    //obsNode->printLocalCoeffs(coeffFile);
    //Node::printAngularSamples(obsLevel);

     auto obsLeaf = dynamic_pointer_cast<Leaf>(obsNode);
     obsLeaf->testMpoleToLocalInLeaf();

    return 0;

}