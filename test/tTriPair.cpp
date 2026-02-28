#include <fstream>
#include "../src/main.h"
#include "../src/clock.h"
#include "../src/config.h"
#include "../src/source.h"
#include "../src/fmm/fmm.h"

using namespace FMM;

extern const Config config("config/config.txt");
extern double k = 0.0;
extern auto states = States();
extern auto t = ClockTimes();

int main() {
    // ===================== Build sources ==================== //
    std::cout << " Importing excitation and sources...\n";

    auto Einc = Exc::importPlaneWaves("config/pwave.txt");
    auto srcs = importSources(Einc);
    size_t nsrcs = srcs.size();

    // ==================== Build nodes ==================== //
    std::cout << " Building FMM tree...\n";

    auto start0 = Clock::now();
    bool isRootLeaf = nsrcs <= config.maxNodeSrcs;
    auto root = std::make_shared<Node>(srcs, 0, nullptr, isRootLeaf);

    root->buildLists();

    std::cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    std::cout << "   # Leaves: " << leaves.size() << '\n';
    std::cout << "   Max node level: " << Node::getMaxLvl() << "\n\n";

    // return 0;

    // ==================== Build nearfield ===================== //
    std::cout << " Building nearfield matrix...\n";

    auto start = Clock::now();
    auto nf = std::make_shared<Nearfield>();
    Time duration_ms = Clock::now() - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // Print tri pairs
    std::cout << Mesh::glTriPairs.size() << " nearfield triangle pairs\n";
    for (const auto& triPair : Mesh::glTriPairs)
        std::cout << "Tri pair: (" 
            << triPair.first.first << ", " 
            << triPair.first.second << ") has "
            << triPair.second.getNumCommonVerts() << " common vertices\n";

    return 0;

}