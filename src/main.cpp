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
    auto root = std::make_shared<Node>(srcs, 0, nullptr, 
        nsrcs <= config.maxNodeSrcs);

    root->buildLists();

    std::cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    std::cout << "   # Leaves: " << leaves.size() << '\n';
    std::cout << "   Max node level: " << Node::getMaxLvl() << "\n\n";

    // ==================== Build nearfield ===================== //
    std::cout << " Building nearfield matrix...\n";

    auto start = Clock::now();
    auto nf = std::make_shared<Nearfield>();
    Time duration_ms = Clock::now() - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build FMM operators =============== //
    std::cout << " Building FMM operators...\n";

    start = Clock::now();
    buildTables();
    root->resizeCoeffs(); // TODO: Hide this call
    duration_ms = Clock::now() - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build radpats ===================== //
    std::cout << " Building radiation patterns...\n";

    start = Clock::now();
    Node::buildRadPats();
    duration_ms = Clock::now() - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Solve iterative FMM ================ //
    std::cout << " Solving w/ FMM...\n";

    constexpr int MAX_ITER = 500;
    constexpr double EPS = 1.0E-6;

    auto solver = std::make_unique<Solver>(srcs, root, nf, MAX_ITER, EPS);

    start = Clock::now();
    solver->solve("curr.txt");
    duration_ms = Clock::now() - start;
    Time duration_ms0 = Clock::now() - start0;
    std::cout << "   FMM total elapsed time: " << duration_ms0.count() << " ms\n\n";

    root->printScatteredField("ff_n"+to_string(nsrcs)+"_selfint.txt", 200);

    if (config.mode == Mode::FMM) return 0;

    // ================== Solve iterative direct ================ //
    states = States(nsrcs);
    resetLeaves();

    start0 = Clock::now();
    root = std::make_shared<Node>(srcs, 0, nullptr, 1);
    root->buildLists();

    std::cout << "\n Building nearfield matrix...\n";

    start = Clock::now();
    nf = std::make_shared<Nearfield>();
    duration_ms = Clock::now() - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    std::cout << " Solving w/ direct...\n";
    solver = std::make_unique<Solver>(srcs, root, nf, MAX_ITER, EPS);

    solver->solve("currDir.txt");
    duration_ms0 = Clock::now() - start0;
    std::cout << "   Direct total elapsed time: " << duration_ms0.count() << " ms\n\n";

    root->printScatteredField("ffDir_n"+to_string(nsrcs)+"_selfint.txt", 200);

    return 0;
}