#include <fstream>
#include "../src/MLFMA.h"
#include "../src/clock.h"
#include "../src/config.h"
#include "../src/fileio.h"

using namespace FMM;

extern auto t = ClockTimes();

int main() {
    // ===================== Read config ==================== //
    std::cout << " Building sources...\n";

    Config config("config/config.txt");

    auto start = Clock::now();
    auto [srcs, Einc] = importFromConfig(config);
    auto end = Clock::now();
    Time duration_ms = end - start;

    auto nsrcs = srcs.size();
    initGlobal(config, Einc, nsrcs);

    // ==================== Test RWG funcs ===================== //
    const auto tri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[0];
    const vec3d obs = { 0.0, 0.0, 5.0 };

    std::cout << std::setprecision(15);

    const auto [scaRad, vecRad, obsProj] = tri.getNearIntegratedRads(obs, 1);
    std::cout << "Numeric :\n   scaRad = " << scaRad << "\n   vecRad = " << vecRad.transpose() << "\n\n";

    const auto [scaRadAnl, vecRadAnl, obsProjAnl] = tri.getNearIntegratedRads(obs, 0);
    std::cout << "Analytic :\n   scaRad = " << scaRadAnl << "\n   vecRad = " << vecRadAnl.transpose() << "\n";

    std::cout << scaRadAnl/scaRad << ' ' << vecRadAnl[0]/vecRad[0] << ' ' << vecRadAnl[1]/vecRad[1] << ' ' << vecRadAnl[2]/vecRad[2] << '\n';

    return 0;
}