#include <fstream>
#include "main.h"
#include "clock.h"
#include "config.h"
#include "source.h"
#include "fmm/fmm.h"

using namespace FMM;

extern const Config config("config/config.txt");
extern auto t = ClockTimes();

int main() {
    // ==================== Import sources ====================== //
    auto Einc = Exc::importPlaneWaves("config/pwave.txt");
    auto srcs = Mesh::importMesh(
        "config/rwg/sph_r5.0_n"+to_string(config.nsrcs), Einc);
    size_t nsrcs = srcs.size();

    // ==================== Build nodes ========================= //
    auto start0 = Clock::now();

    bool isRootLeaf = nsrcs <= config.maxNodeSrcs;
    auto root = std::make_shared<Node>(srcs, 0, nullptr, isRootLeaf);

    root->buildLists();

    // ==================== Build nearfield ===================== //
    auto nf = std::make_shared<Nearfield>();

    // ==================== Build FMM operators ================= //
    if (!isRootLeaf) buildTables();

    // ==================== Build expansions ==================== //
    root->resizeCoeffs(); // TODO: Hide this call
    if(!isRootLeaf) buildRadPats();

    // ==================== Solve iterative FMM ================= //
    const int MAX_ITER = nsrcs;
    constexpr double EPS = 1.0E-6;

    auto solver = std::make_unique<GMRES>(srcs, nf, root, EPS, MAX_ITER);
    solver->solve("sol.txt");

    Time duration_ms0 = Clock::now() - start0;
    std::cout << " FMM total elapsed time: " << duration_ms0.count() << " ms\n\n";

    Mesh::printScattered(srcs, "ff_n"+to_string(nsrcs)+".txt", 200);

    if (config.mode == Mode::FMM) return 0;

    // ================== Solve iterative direct ================ //
    start0 = Clock::now();

    resetNodes();
    root = std::make_shared<Node>(srcs, 0, nullptr, 1);
    root->buildLists();

    nf = std::make_shared<Nearfield>();

    solver = std::make_unique<GMRES>(srcs, nf, root, EPS, MAX_ITER);
    solver->solve("solDir.txt");
    //auto solverDir = std::make_unique<Direct>(srcs, nf);
    //solverDir->solve("solDir.txt");

    duration_ms0 = Clock::now() - start0;
    std::cout << " Direct total elapsed time: " << duration_ms0.count() << " ms\n\n";

    Mesh::printScattered(srcs, "ffDir_n"+to_string(nsrcs)+".txt", 200);

    return 0;
}