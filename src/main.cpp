#include <fstream>
#include "main.h"
#include "clock.h"
#include "config.h"
#include "source.h"
#include "states.h"
#include "fmm/fmm.h"

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

    states = States(nsrcs);

    // ==================== Build nodes ==================== //
    std::cout << " Building FMM tree...\n";

    auto start0 = Clock::now();
    bool isRootLeaf = nsrcs <= config.maxNodeSrcs;
    auto root = std::make_shared<Node>(srcs, 0, nullptr, isRootLeaf);

    root->buildLists();

    std::cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    std::cout << "   # Leaves: " << leaves.size() << '\n';
    std::cout << "   Max node level: " << Node::getMaxLvl() << "\n\n";

    // ==================== Build nearfield ===================== //
    std::cout << " Building nearfield matrix...     ";

    auto start = Clock::now();
    auto nf = std::make_shared<Nearfield>();
    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";

    // ==================== Build FMM operators =============== //
    std::cout << " Building FMM operators...        ";

    start = Clock::now();
    if (!isRootLeaf) buildTables();
    duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";

    // ==================== Build expansions ==================== //
    std::cout << " Building plane wave expansions...";

    start = Clock::now();
    root->resizeCoeffs(); // TODO: Hide this call
    if(!isRootLeaf) buildRadPats();
    duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";

    // ==================== Solve iterative FMM ================ //
    std::cout << " Solving with FMM...              ";

    constexpr int MAX_ITER = 500;
    constexpr double EPS = 1.0E-6;

    auto solver = std::make_unique<Solver>(srcs, root, nf, MAX_ITER, EPS);
    solver->solve("curr.txt");

    Time duration_ms0 = Clock::now() - start0;
    std::cout << " FMM total elapsed time: " << duration_ms0.count() << " ms\n\n";

    Mesh::printScattered(srcs, "ff_n"+to_string(nsrcs)+".txt", 200);

    if (config.mode == Mode::FMM) return 0;

    // ================== Solve iterative direct ================ //
    states = States(nsrcs);
    resetLeaves();

    start0 = Clock::now();
    root = std::make_shared<Node>(srcs, 0, nullptr, 1);
    root->buildLists();

    std::cout << "\n Building nearfield matrix...     ";

    start = Clock::now();
    nf = std::make_shared<Nearfield>();
    duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";

    std::cout << " Solving with direct...           ";
    solver = std::make_unique<Solver>(srcs, root, nf, MAX_ITER, EPS);
    solver->solve("currDir.txt");

    duration_ms0 = Clock::now() - start0;
    std::cout << " Direct total elapsed time: " << duration_ms0.count() << " ms\n\n";

    Mesh::printScattered(srcs, "ffDir_n"+to_string(nsrcs)+".txt", 200);

    return 0;
}