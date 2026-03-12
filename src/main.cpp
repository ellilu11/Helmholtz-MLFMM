#include <fstream>
#include "main.h"
#include "clock.h"
#include "config.h"
#include "dipole.h"
#include "fmm/fmm.h"

extern const Config config("config/config.txt");

void mainLoop(const SrcVec& srcs, bool doFMM, bool doIter = true) {
    size_t nsrcs = srcs.size();
    std::string method = doFMM ? "FMM" : "Direct";

    // ==================== Build nodes ========================= //
    auto start = Clock::now();

    bool isRootLeaf = nsrcs <= config.maxNodeSrcs || !doFMM;
    auto root = std::make_shared<FMM::Node>(srcs, 0, nullptr, isRootLeaf);
    root->buildLists();
   
    // ==================== Build nearfield ===================== //
    auto nf = std::make_unique<FMM::Nearfield>(nsrcs);

    // ==================== Build FMM quantities ================ //
    if (!isRootLeaf) {
        FMM::buildLevels();
        root->resizeCoeffs(); // TODO: Hide this call
        FMM::buildRadPats();
    }

    // ==================== Solve for current =================== //
    std::unique_ptr<Solver> solver;
    if (doFMM || (!doFMM && doIter))
        solver = std::make_unique<GMRES>(srcs, std::move(nf), root, 1.0E-6, nsrcs);
    else 
        solver = std::make_unique<Direct>(srcs, std::move(nf));

    solver->solve(doFMM ? "sol.txt" : "solDir.txt");

    Time duration_ms = Clock::now() - start;
    std::cout << " " + method + " total elapsed time : " << duration_ms.count() << " ms\n\n";

    // ==================== Compute scattered field ============= //
    Mesh::printScattered(srcs,
        "out/ff/px_k1.0z_r5.0_efie",
        (doFMM ? "ff_n" : "ffDir_n")+std::to_string(nsrcs)+".txt", 200);
    //Mesh::printScattered(srcs,
    //    "out/ff/py_km1.0x_plate",
    //    std::string(doFMM ? "ff_g" : "ffDir_g")+"3.00.txt", 100);
}

int main() {
    auto Einc = Exc::importPlaneWaves("config/pwave.txt");
    auto srcs = Mesh::importMesh(
        "config/rwg/sph_r5.0/sph_r5.0_n"+std::to_string(config.nsrcs), Einc);
    // Mesh::printNormals("out/nhats.txt");

    constexpr bool doIter = true;
    switch (config.mode) {
        case Mode::FMM: mainLoop(srcs, true); break;
        case Mode::DIR: mainLoop(srcs, false, doIter); break;
        case Mode::FMMDIR: {
            mainLoop(srcs, true);
            FMM::reset();
            mainLoop(srcs, false, doIter);
        } break;
    }

    return 0;
}