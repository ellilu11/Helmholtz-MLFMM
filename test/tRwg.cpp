#include <fstream>
#include "../src/main.h"
#include "../src/clock.h"
#include "../src/config.h"
#include "../src/fileio.h"
#include "../src/mesh/triangle.h"

using namespace FMM;

extern const Config config("config/config.txt");
extern double k = 0.0;
extern auto t = ClockTimes();

void testSingularityExtraction(const Mesh::Triangle& tri, const vec3d& obs) {
    std::cout << std::setprecision(15);

    const auto [scaRad, vecRad] = tri.getNearIntegrated(obs, 1);
    std::cout << "Numeric :\n   scaRad = " << scaRad << "\n   vecRad = " << vecRad.transpose() << "\n\n";

    const auto [scaRadAnl, vecRadAnl] = tri.getNearIntegrated(obs, 0);
    std::cout << "Analytic :\n   scaRad = " << scaRadAnl << "\n   vecRad = " << vecRadAnl.transpose() << "\n";

    std::cout << scaRadAnl/scaRad << ' ' 
        << vecRadAnl[0]/vecRad[0] << ' ' 
        << vecRadAnl[1]/vecRad[1] << ' ' 
        << vecRadAnl[2]/vecRad[2] << '\n';
}

int main() {
    // ===================== Read config ==================== //
    std::cout << " Building sources...\n";

    auto start = Clock::now();
    auto [srcs, Einc] = importFromConfig(config);
    auto end = Clock::now();
    Time duration_ms = end - start;

    auto nsrcs = srcs.size();
    initGlobal(Einc, nsrcs);

    // Singularity extraction test
    const auto tri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[0];
    const vec3d obs = { 0.0, 0.0, 5.0 };
    testSingularityExtraction(tri, obs);
    //

    /* Nearfield integration test
    auto rwg0 = dynamic_pointer_cast<Mesh::RWG>(srcs[0]);
    auto rwg1 = dynamic_pointer_cast<Mesh::RWG>(srcs[1]);

    auto selfRad = rwg0->getIntegratedRad(rwg0);
    auto pairRad = rwg1->getIntegratedRad(rwg0);

    std::cout << std::setprecision(15)
        << "Self intRad: " << selfRad << '\n'
        << "Pair intRad: " << pairRad << '\n'
        << "Diff: " << pairRad - selfRad << '\n';
    */

    return 0;
}

