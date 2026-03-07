#include <fstream>
#include "../src/main.h"
#include "../src/clock.h"
#include "../src/config.h"
#include "../src/source.h"
#include "../src/mesh/triangle.h"

using namespace FMM;

extern const Config config("config/config.txt");
extern auto t = ClockTimes();

void testInvRNumVsAnl(const Mesh::Triangle& tri, const vec3d& obs) {
    std::cout << std::setprecision(15);

    const auto [scaRad, vecRad] = tri.getIntegratedInvRcubed(obs, 1);
    std::cout << "Numeric :\n   scaRad = " << scaRad << "\n   vecRad = " << vecRad.transpose() << "\n\n";

    const auto [scaRadAnl, vecRadAnl] = tri.getIntegratedInvRcubed(obs, 0);
    std::cout << "Analytic :\n   scaRad = " << scaRadAnl << "\n   vecRad = " << vecRadAnl.transpose() << "\n";

    std::cout << scaRadAnl/scaRad << ' ' 
        << vecRadAnl[0]/vecRad[0] << ' ' 
        << vecRadAnl[1]/vecRad[1] << ' ' 
        << vecRadAnl[2]/vecRad[2] << '\n';
}

void testSingularEFIE(const Mesh::Triangle& obsTri, const Mesh::Triangle& srcTri) {
    assert(config.alpha == 1.0); // only EFIE part should be tested here

    double k2 = config.k * config.k;
    vec3d vobs = obsTri.getVerts()[0];
    vec3d vsrc = srcTri.getVerts()[0];

    pair2i pair = std::minmax(obsTri.getIdx(), srcTri.getIdx());
    Mesh::TriPair triPair(pair);

    /* Double numeric integration of 1/R 
    cmplx radNumNum = 0.0;
    for (const auto& [obs, obsWeight] : obsTri.getQuads()) {
        for (const auto& [src, srcWeight] : srcTri.getQuads()) {
            double r = (obs-src).norm();
            assert(!Math::fzero(r));

            cmplx G = 1.0 / r;
            radNumNum += ((obs-vobs).dot(src-vsrc) - 4.0/k2) * G
                * obsWeight * srcWeight;
        }
    }*/

    // Double numeric integration of 1/R using precomputed moments
    const auto& [m00, m10, m01, m11] = triPair.momentsEFIE;
    vec3d v0 = (obsTri.getIdx() <= srcTri.getIdx()) ? vobs : vsrc;
    vec3d v1 = (obsTri.getIdx() <= srcTri.getIdx()) ? vsrc : vobs;
    cmplx radNumNum = m11 - v1.dot(m10) - v0.dot(m01) + (v0.dot(v1) - 4.0/k2)*m00;
    //

    // Numeric-analytic integration of 1/R
    cmplx radNumAnl = (
        obsTri.getSingularEFIE(srcTri, triPair, vobs, vsrc) +
        srcTri.getSingularEFIE(obsTri, triPair, vsrc, vobs)) / 2.0;

    std::cout << std::setprecision(15)
        << "NumNum EFIE: " << radNumNum << '\n'
        << "NumAnl EFIE: " << radNumAnl << '\n'
        << "Diff: " << radNumNum - radNumAnl << '\n';
}


void testSingularMFIE(const Mesh::Triangle& obsTri, const Mesh::Triangle& srcTri) {
    assert(config.alpha == 0.0); // only MFIE part should be tested here

    double k2 = config.k * config.k;
    vec3d vobs = obsTri.getVerts()[0];
    vec3d vsrc = srcTri.getVerts()[0];

    pair2i pair = std::minmax(obsTri.getIdx(), srcTri.getIdx());
    Mesh::TriPair triPair(pair);

    /* Double numeric integration of k^2/2 * 1/R + 1/R^3
    cmplx radNumNum = 0.0;
    for (const auto& [obs, obsWeight] : obsTri.getQuads()) {
        for (const auto& [src, srcWeight] : srcTri.getQuads()) {
            vec3d R = obs-src;
            double r = R.norm(), r2 = r*r, r3 = r*r2;
            assert(!Math::fzero(r));

            vec3cd gradG = R/r3 * (0.5*k2*r2 + 1.0);

            // minus sign from flipping J x gradG to gradG x J
            radNumNum -= conj((obs-vobs).dot(gradG.cross(src-vsrc))) // Hermitian cross!
                * obsWeight * srcWeight;
        }
    }
    */

    // Double numeric integration of k^2/2 * 1/R + 1/R^3 using precomputed moments
    const auto& [m00, m10, m01, m11] = triPair.momentsMFIE;
    vec3d v0 = (obsTri.getIdx() <= srcTri.getIdx()) ? vobs : vsrc;
    vec3d v1 = (obsTri.getIdx() <= srcTri.getIdx()) ? vsrc : vobs;
    cmplx radNumNum = m11 - v1.dot(m10) - v0.dot(m01) + (v1.cross(v0)).dot(m00);
    //

    // Numeric-analytic integration of k^2/2 * 1/R + 1/R^3
    cmplx radNumAnl = -obsTri.getSingularMFIE(srcTri, triPair, vobs, vsrc);

    std::cout << std::setprecision(15)
        << "NumNum MFIE: " << radNumNum << '\n'
        << "NumAnl MFIE: " << radNumAnl << '\n'
        << "Diff: " << radNumNum - radNumAnl << '\n';
}

int main() {
    // ===================== Build sources ==================== //
    std::cout << " Building sources...\n";

    const auto Einc = Exc::importPlaneWaves("config/pwave.txt");
    const auto srcs = importSources(Einc);
    size_t nsrcs = srcs.size();

    /* Numeric vs analytic 1/R or 1/R^3 integration test
    const auto tri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[0];
    const vec3d obs = { -1.00000, 2.0, 0.0 };
    testInvRNumVsAnl(tri, obs);
    */

    // Double numeric vs numeric-analytic singular EFIE or MFIE test
    auto srcTri = dynamic_pointer_cast<Mesh::RWG>(srcs[500])->getTris()[1];
    auto obsTri = dynamic_pointer_cast<Mesh::RWG>(srcs[5])->getTris()[0];
    testSingularMFIE(obsTri, srcTri);
    std::cout << " (obsTri.idx, srcTri.idx) = (" << obsTri.getIdx() << ", " << srcTri.getIdx() << ")\n";
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

