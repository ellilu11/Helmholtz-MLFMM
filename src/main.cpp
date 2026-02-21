#include <fstream>
#include "main.h"
#include "clock.h"
#include "config.h"
#include "fileio.h"
#include "fmm/fmm.h"

using namespace FMM;

extern const Config config("config/config.txt");
extern double k = 0.0;
extern auto t = ClockTimes();

int main() {
    // ===================== Build sources ==================== //
    std::cout << " Building sources...\n";

    auto [srcs, Einc] = importFromConfig(config);
    auto nsrcs = srcs.size();

    initGlobal(Einc, nsrcs);

    // ==================== Build nodes ==================== //
    std::cout << " Building nodes...\n";

    shared_ptr<Node> root;
    if (nsrcs > config.maxNodeSrcs)
        root = std::make_shared<Stem>(srcs, 0, nullptr);
    else
        root = std::make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    std::cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    std::cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    std::cout << "   Max node level: " << Node::getMaxLvl() << "\n\n";

    // ==================== Build nearfield ===================== //
    std::cout << " Building nearfield interactions...\n";

    auto start = Clock::now();
    Leaf::buildNearRads();
    auto end = Clock::now();
    Time duration_ms = end - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build FMM operators =============== //
    std::cout << " Building FMM operators...\n";

    start = Clock::now();
    buildTables();
    root->resizeCoeffs(); // TODO: Hide this call
    end = Clock::now();
    duration_ms = end - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Build radpats ===================== //
    std::cout << " Building radiation patterns...\n";

    start = Clock::now();
    Leaf::buildRadPats();
    end = Clock::now();
    duration_ms = end - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    // ==================== Solve iterative FMM ================ //
    std::cout << " Solving w/ FMM...\n";

    constexpr int MAX_ITER = 500;
    constexpr double EPS = 1.0E-6;

    auto solver = std::make_unique<Solver>(srcs, root, MAX_ITER, EPS,
        lvec, rvec, currents);

    start = Clock::now();
    solver->solve("curr.txt");
    end = Clock::now();
    duration_ms = end - start;
    std::cout << "   Total elapsed time: " << duration_ms.count() << " ms\n\n";

    // root->printFarFld("ff.txt");

    if (!config.evalDirect) return 0;

    // ================== Solve iterative direct ================ //
    initGlobal(Einc, nsrcs);

    root = std::make_shared<Leaf>(srcs, 0, nullptr);
    root->buildLists();

    std::cout << " Building nearfield interactions...\n";

    start = Clock::now();
    Leaf::buildNearRads();
    end = Clock::now();

    duration_ms = end - start;
    std::cout << "   Elapsed time: " << duration_ms.count() << " ms\n\n";

    std::cout << " Solving w/ direct...\n";
    solver = std::make_unique<Solver>(srcs, root, MAX_ITER, EPS,
        lvec, rvec, currents);

    solver->solve("currDir.txt");

    return 0;
}