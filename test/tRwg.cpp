#include <fstream>
#include "../src/main.h"
#include "../src/clock.h"
#include "../src/config.h"
#include "../src/source.h"
#include "../src/mesh/triangle.h"

using namespace FMM;

extern const Config config("config/config.txt");

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

void testMomentsEFIE(const Mesh::Triangle& obsTri, const Mesh::Triangle& srcTri) {
    assert(config.alpha == 1.0); // only EFIE part should be tested here

    double k = config.k, k2 = k*k;
    vec3d vobs = obsTri.getVerts()[0];
    vec3d vsrc = srcTri.getVerts()[0];

    pair2i pair = std::minmax(obsTri.getIdx(), srcTri.getIdx());
    Mesh::TriPair triPair(pair);

    cmplx rad = 0.0;
    for (const auto& [obs, obsWeight] : obsTri.getQuads()) {
        for (const auto& [src, srcWeight] : srcTri.getQuads()) {
            double r = (obs-src).norm();
            assert(!Math::fzero(r));

            cmplx G = obsWeight*srcWeight / (4.0*PI);
            if (triPair.nCommon >= 2)
                G *= (Math::fzero(r) ? iu*k : (exp(iu*k*r)-1.0) / r);
            else {
                assert(!Math::fzero(r));
                G *= exp(iu*k*r) / r;
            }

            rad += ((obs-vobs).dot(src-vsrc) - 4.0/k2) * G;
        }
    }

    const auto& [m00, m10, m01, m11] = triPair.momentsEFIE;
    vec3d v0 = (obsTri.getIdx() <= srcTri.getIdx()) ? vobs : vsrc;
    vec3d v1 = (obsTri.getIdx() <= srcTri.getIdx()) ? vsrc : vobs;
    cmplx radMoments = m11 - v1.dot(m10) - v0.dot(m01) + (v0.dot(v1) - 4.0/k2)*m00;

    std::cout << std::setprecision(15)
        << "Direct: " << rad << '\n'
        << "Moments: " << radMoments << '\n'
        << "Diff: " << rad - radMoments << '\n';
}

void testSingularEFIE(const Mesh::Triangle& obsTri, const Mesh::Triangle& srcTri, std::ofstream& os) {
    assert(config.alpha == 1.0); // only EFIE part should be tested here
    if (obsTri.getIdx() == srcTri.getIdx()) return;

    double k2 = config.k * config.k;
    vec3d vobs = obsTri.getVerts()[0];
    vec3d vsrc = srcTri.getVerts()[0];

    pair2i pair = std::minmax(obsTri.getIdx(), srcTri.getIdx());
    Mesh::TriPair triPair(pair);

    // Double numeric integration of 1/R 
    double radNumNum = 0.0;
    for (const auto& [obs, obsWeight] : obsTri.getQuads()) {
        for (const auto& [src, srcWeight] : srcTri.getQuads()) {
            double r = (obs-src).norm();

            double G = 1.0 / (4.0*PI*r);
            radNumNum += ((obs-vobs).dot(src-vsrc) - 4.0/k2) * G
                * obsWeight * srcWeight;
        }
    }

    // Numeric-analytic integration of 1/R
    double radNumAnl = (
        obsTri.getSingularEFIE(srcTri, triPair, vobs, vsrc) +
        srcTri.getSingularEFIE(obsTri, triPair, vsrc, vobs)) / 2.0;

    //std::cout << std::setprecision(15)
    //    << "NumNum EFIE: " << radNumNum << ' '
    //    << "NumAnl EFIE: " << radNumAnl << ' '
    //    << "Diff: " << radNumNum - radNumAnl << '\n';

    //std::cout << std::setprecision(9) << " # common verts: " << triPair.nCommon << ' '
    //    << "Rel diff: " << (radNumNum - radNumAnl) / radNumAnl << '\n';

    os << std::setprecision(9)
        << (obsTri.getCenter() - srcTri.getCenter()).norm() << ' '
        << (radNumNum - radNumAnl) / radNumAnl << '\n';
}

void testSelfIntegrated(const Mesh::Triangle& obsTri, const Mesh::Triangle& srcTri) {
    assert(config.alpha == 1.0); // only EFIE part should be tested here

    std::cout << std::setprecision(9);
    std::cout << "obsTri #" << obsTri.getIdx() << '\n';
    std::cout << "srcTri #" << srcTri.getIdx() << '\n';

    vec3d vobs = obsTri.getVerts()[0];
    vec3d vsrc = obsTri.getVerts()[1];

    pair2i pair = std::minmax(obsTri.getIdx(), srcTri.getIdx());
    Mesh::TriPair triPair(pair);

    double selfRad = obsTri.getDoubleSelfIntegratedInvR(vobs, vsrc);
    double pairRad = obsTri.getSingularEFIE(srcTri, triPair, vobs, vsrc);
        //(obsTri.getSingularEFIE(srcTri, triPair, vobs, vsrc) +
        //    srcTri.getSingularEFIE(obsTri, triPair, vsrc, vobs)) / 2.0;

    std::cout << "Z diff: " <<
        srcTri.getVerts()[0][2] - obsTri.getVerts()[0][2] << '\n';

    std::cout << std::setprecision(15)
        << "Self intRad: " << selfRad << '\n'
        << "Pair intRad: " << pairRad << '\n'
        << "Diff: " << pairRad - selfRad << '\n';
}

int main() {
    // ==================== Import sources ====================== //
    auto Einc = Exct::importPlaneWaves("config/pwave.txt");
    auto srcs = Mesh::importMesh(
        "config/rwg/test/n"+std::to_string(config.nsrcs)+"self", Einc);
    size_t nsrcs = srcs.size();

    /* ==================== Build nodes ========================= //
    auto start = Clock::now();

    bool isRootLeaf = nsrcs <= config.maxNodeSrcs;
    auto root = std::make_shared<FMM::Node>(srcs, 0, nullptr, isRootLeaf);
    root->buildLists();

    // ==================== Build nearfield ===================== //
    auto nf = std::make_unique<FMM::Nearfield>(nsrcs);
    */

    // ===================== Test ========================== //
    /*
    std::ofstream file("out/test/mesh/singular_efie_q13.txt");

    for (const auto& triMap : Mesh::glTriPairs) {
        const auto& [obsTri, srcTri] = triMap.second.getTriPair();

        testSingularEFIE(obsTri, srcTri, file);
        testSingularEFIE(srcTri, obsTri, file);
    }*/

    /* Numeric vs analytic 1/R or 1/R^3 integration test
    const auto tri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[0];
    const vec3d obs = { -1.00000, 2.0, 0.0 };
    testInvRNumVsAnl(tri, obs);
    */

    /* Double numeric vs numeric-analytic singular EFIE or MFIE test
    auto srcTri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[0];
    auto obsTri = dynamic_pointer_cast<Mesh::RWG>(srcs[1])->getTris()[1];
    testSingularMFIE(obsTri, srcTri);
    std::cout << "(obsTri, srcTri) = (" << obsTri.getIdx() << ", " << srcTri.getIdx() << ")\n";
    */

    /* Mass matrix numeric vs analytic test
    auto tri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[0];
    testMassMatrix(tri);
    */

    // Near/self integration test
    auto obsTri = dynamic_pointer_cast<Mesh::RWG>(srcs[0])->getTris()[0];
    auto srcTri = dynamic_pointer_cast<Mesh::RWG>(srcs[1])->getTris()[0];
    testSelfIntegrated(obsTri, srcTri);

    return 0;
}

