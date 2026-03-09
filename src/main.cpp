#include <fstream>
#include "main.h"
#include "clock.h"
#include "config.h"
#include "source.h"
#include "fmm/fmm.h"

using namespace FMM;

extern const Config config("config/config.txt");
extern auto t = ClockTimes();

void mainLoop(const SrcVec& srcs, bool doFMM, bool doGMRES = true) {
    size_t nsrcs = srcs.size();
    std::string method = doFMM ? "FMM" : "Direct";

    std::string ieStr = getIEStr(config.ie);
    std::transform(ieStr.begin(), ieStr.end(), ieStr.begin(), ::tolower);

    // ==================== Build nodes ========================= //
    auto start = Clock::now();

    bool isRootLeaf = nsrcs <= config.maxNodeSrcs || !doFMM;
    auto root = std::make_shared<Node>(srcs, 0, nullptr, isRootLeaf);

    root->buildLists();

    // ==================== Build nearfield ===================== //
    auto nf = std::make_shared<Nearfield>();

    // ==================== Build FMM operators ================= //
    if (!isRootLeaf) buildTables();

    // ==================== Build expansions ==================== //
    root->resizeCoeffs(); // TODO: Hide this call
    if (!isRootLeaf) buildRadPats();

    // ==================== Solve for current =================== //
    std::unique_ptr<Solver> solver;
    if (doFMM || (!doFMM && doGMRES))
        solver = std::make_unique<GMRES>(srcs, nf, root, 1.0E-6, config.maxIter);
    else 
        solver = std::make_unique<Direct>(srcs, nf);

    solver->solve(doFMM ? "sol.txt" : "solDir.txt");

    Time duration_ms = Clock::now() - start;
    std::cout << " " + method + " total elapsed time : " << duration_ms.count() << " ms\n\n";

    // ==================== Compute scattered field ============= //
    Mesh::printScattered(srcs,
        (doFMM ? "ff_n" : "ffDir_n")+to_string(nsrcs)+"_"+ieStr+".txt", 200);
}

int main() {
    auto Einc = Exc::importPlaneWaves("config/pwave.txt");
    auto srcs = Mesh::importMesh(
        "config/rwg/sph_r5.0_n"+to_string(config.nsrcs), Einc);

    constexpr bool doGMRES = false;
    switch (config.mode) {
        case Mode::FMM: mainLoop(srcs, true); break;
        case Mode::DIR: mainLoop(srcs, false, doGMRES); break;
        case Mode::FMMDIR: {
            mainLoop(srcs, true);
            resetNodes();
            mainLoop(srcs, false, doGMRES);
        }
        break;
    }

    return 0;
}