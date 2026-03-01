#include <fstream>
#include "../src/main.h"
#include "../src/clock.h"
#include "../src/config.h"
#include "../src/source.h"
#include "../src/mesh/triangle.h"

using namespace FMM;

extern const Config config("config/config.txt");
extern double k = 0.0;
extern auto states = States();
extern auto t = ClockTimes();

void testInvRNumVsAnl(const Mesh::Triangle& tri, const vec3d& obs) {
    std::cout << std::setprecision(15);

    const auto [scaRad, vecRad] = tri.getIntegratedInvR(obs, 1);
    std::cout << "Numeric :\n   scaRad = " << scaRad << "\n   vecRad = " << vecRad.transpose() << "\n\n";

    const auto [scaRadAnl, vecRadAnl] = tri.getIntegratedInvR(obs, 0);
    std::cout << "Analytic :\n   scaRad = " << scaRadAnl << "\n   vecRad = " << vecRadAnl.transpose() << "\n";

    std::cout << scaRadAnl/scaRad << ' ' 
        << vecRadAnl[0]/vecRad[0] << ' ' 
        << vecRadAnl[1]/vecRad[1] << ' ' 
        << vecRadAnl[2]/vecRad[2] << '\n';
}

void testInvR(const Mesh::Triangle& obsTri, const Mesh::Triangle& srcTri) {
    
    const auto& vobs = obsTri.getVerts()[0];
    const auto& vsrc = srcTri.getVerts()[0];
    const auto& vsrcproj = srcTri.proj(vsrc);
    
    const int nCommon = obsTri.getNumCommonVerts(srcTri);
    std::cout << "nCommon = " << nCommon << '\n';
    
    cmplx nearRad = 0.0, fullRad = 0.0;

    /*
    for (const auto& [obs, obsWeight] : obsTri.getQuads()) {

        for (const auto& [src, srcWeight] : srcTri.getQuads()) {
            const double r = (obs-src).norm();
            cmplx G = 1.0 / r;

            fullRad += ((obs-vobs).dot(src-vsrc) - 4.0 / (k*k)) * G
                * obsWeight * srcWeight;
        }

        const auto& obsProj = srcTri.proj(obs);
        const auto& [scaRad, vecRad] = srcTri.getIntegratedInvR(obs);
        nearRad +=
            ((obs-vobs).dot(vecRad+(obsProj-vsrcproj)*scaRad) - 4.0/(k*k)*scaRad)
            * obsWeight;
    }*/

    std::cout << std::setprecision(15)
        << "Full Rad: " << fullRad << '\n'
        << "Near Rad: " << nearRad << '\n'
        << "Diff: " << fullRad-nearRad << '\n';
}

int main() {
    // ===================== Build sources ==================== //
    std::cout << " Building sources...\n";

    const auto Einc = Exc::importPlaneWaves("config/pwave.txt");
    const auto srcs = importSources(Einc);
    size_t nsrcs = srcs.size();

    /* Numeric vs analytic 1/R near integration test
    const auto tri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[0];
    const vec3d obs = { -1.00000, 2.0, 0.0 };
    testInvRNumVsAnl(tri, obs);
    */

    // Numeric full vs analytic near 1/R integration test
    const auto obsTri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[1];
    const auto srcTri = dynamic_pointer_cast<Mesh::RWG>(srcs[1])->getTris()[0];
    testInvR(obsTri, srcTri);
    //

    /* Near/self integration test
    auto rwg0 = dynamic_pointer_cast<Mesh::RWG>(srcs[0]);
    auto rwg1 = dynamic_pointer_cast<Mesh::RWG>(srcs[1]);

    auto selfRad = rwg0->getIntegratedRad(rwg0);
    auto pairRad = rwg1->getIntegratedRad(rwg0);

    //std::cout << std::setprecision(15)
    //    << "Self intRad: " << selfRad << '\n'
    //    << "Pair intRad: " << pairRad << '\n'
    //    << "Diff: " << pairRad - selfRad << '\n';
    
    std::cout << std::setprecision(15)
        << rwg1->getTris()[0].getVerts()[0][2] - rwg0->getTris()[0].getVerts()[0][2] << ' '
        << (pairRad-selfRad).real() << ' ' 
        << pairRad.real() << '\n';
    */

    return 0;
}

