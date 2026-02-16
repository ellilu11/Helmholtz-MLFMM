#include <fstream>
#include "../src/MLFMA.h"
#include "../src/clock.h"
#include "../src/config.h"
#include "../src/fileio.h"
#include "../src/sources/mesh/triangle.h"

using namespace FMM;

extern auto t = ClockTimes();
extern bool doDirFar = false;

/*
void testSingularityExtraction(const Triangle& tri, const vec3d& obs) {

    std::cout << std::setprecision(15);

    const auto [scaRad, vecRad, obsProj] = tri.getNearIntegrated(obs, 1);
    std::cout << "Numeric :\n   scaRad = " << scaRad << "\n   vecRad = " << vecRad.transpose() << "\n\n";

    const auto [scaRadAnl, vecRadAnl, obsProjAnl] = tri.getNearIntegrated(obs, 0);
    std::cout << "Analytic :\n   scaRad = " << scaRadAnl << "\n   vecRad = " << vecRadAnl.transpose() << "\n";

    std::cout << scaRadAnl/scaRad << ' ' 
        << vecRadAnl[0]/vecRad[0] << ' ' 
        << vecRadAnl[1]/vecRad[1] << ' ' 
        << vecRadAnl[2]/vecRad[2] << '\n';
}
*/

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

    /* Singularity extraction test
    const auto tri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[0];
    const vec3d obs = { 0.0, 0.0, 5.0 };
    testSingularityExtraction(tri, obs);
    */

    // Nearfield integration test
    double k = Einc->wavenum;

    auto obs = dynamic_pointer_cast<Mesh::RWG>(srcs[1]);
    auto src = dynamic_pointer_cast<Mesh::RWG>(srcs[0]);

    obs->getIntegratedRad(src);

    /*
    int iObsTri = 1, iSrcTri = 0;

    auto obsTri = obs->getTris()[iObsTri];
    auto srcTri = src->getTris()[iSrcTri];

    vec3d obsNC = obs->getVertsNC()[iObsTri];
    vec3d srcNC = src->getVertsNC()[iSrcTri];
    //std::cout << "obsNC: " << obsNC.transpose() << "\n"
    //          << "srcNC: " << srcNC.transpose() << "\n";

    auto [X0, X1, X2] = srcTri.getVerts();
    std::cout << "srcTri verts:\n   " << X0.transpose() << "\n   " << X1.transpose() << "\n   " << X2.transpose() << "\n";

    auto [Y0, Y1, Y2] = obsTri.getVerts();
    std::cout << "obsTri verts:\n   " << Y0.transpose() << "\n   " << Y1.transpose() << "\n   " << Y2.transpose() << "\n";

    const int nCommon = obsTri.getNumCommonVerts(srcTri);
    assert(!nCommon); // No common vertices

    cmplx intRad = 0.0, intRadFull = 0.0;
    for (const auto& [obs, obsWeight] : obsTri.getQuads()) {
        for (const auto& [src, srcWeight] : srcTri.getQuads()) {
            // Using bilinear EFIE
            double r = (obs-src).norm();
            cmplx G = exp(iu*k*r) / r;
            intRad +=
                ((obs-obsNC).dot(src-srcNC) - 4.0/(k*k)) * G
                * obsWeight * srcWeight
                * Math::sign(iObsTri) * Math::sign(iSrcTri);

            // Using full EFIE
            const auto& dyadic = Math::dyadicG(obs-src, k);
            intRadFull +=
                (obs-obsNC).dot(dyadic*(src-srcNC))
                * obsWeight * srcWeight
                * Math::sign(iObsTri) * Math::sign(iSrcTri);
        }
    }

    std::cout << std::setprecision(15)
        << "intRad (bilinear): " << intRad << "\n"
        << "intRad (full): " << intRadFull << "\n";
    */

    return 0;
}