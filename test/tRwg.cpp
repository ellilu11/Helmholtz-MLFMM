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

            cmplx G;
            if (triPair.nCommon >= 2)
                G = (Math::fzero(r) ? iu*k : (exp(iu*k*r)-1.0) / r);
            else {
                assert(!Math::fzero(r));
                G = exp(iu*k*r) / r;
            }

            rad += ((obs-vobs).dot(src-vsrc) - 4.0/k2) * G
                * obsWeight * srcWeight;
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

void testMomentsMFIE(const Mesh::Triangle& obsTri, const Mesh::Triangle& srcTri) {
    assert(config.alpha == 0.0); // only MFIE part should be tested here

    double k = config.k, k2 = k*k;
    vec3d vobs = obsTri.getVerts()[0];
    vec3d vsrc = srcTri.getVerts()[0];
    vec3d nhat = obsTri.getNormal();

    pair2i pair = std::minmax(obsTri.getIdx(), srcTri.getIdx());
    Mesh::TriPair triPair(pair);

    cmplx rad = 0.0;
    for (const auto& [obs, obsWeight] : obsTri.getQuads()) {
        for (const auto& [src, srcWeight] : srcTri.getQuads()) {
            vec3d R = obs-src;
            double r = R.norm(), r2 = r*r, r3 = r*r2;
            assert(!Math::fzero(r));

            vec3cd gradG = R / r3;
            if (triPair.nCommon == 2)
                gradG = gradG * ((-1.0+iu*k*r)*exp(iu*k*r) + 0.5*k*k*r2 + 1.0);
            else if (triPair.nCommon < 2)
                gradG = gradG * (-1.0+iu*k*r)*exp(iu*k*r);

            rad -= (obs-vobs).dot(nhat.cross(gradG.cross(src-vsrc)))
                * obsWeight * srcWeight;
        }
    }

    const auto& [m000, m001, m10, m01, m11] = 
        (obsTri.getIdx() <= srcTri.getIdx()) ? triPair.momentsMFIE : triPair.momentsMFIE2;
    cmplx radMoments = 
        m11 - vsrc.dot(m10) - vobs.dot(m01) + (vobs.dot(vsrc))*m000 + nhat.dot(vsrc)*vobs.dot(m001);

    std::cout << std::setprecision(15)
        << "Direct: " << rad << '\n'
        << "Moments: " << radMoments << '\n'
        << "Diff: " << rad - radMoments << '\n';
}

void testSingularEFIE(const Mesh::Triangle& obsTri, const Mesh::Triangle& srcTri) {
    assert(config.alpha == 1.0); // only EFIE part should be tested here

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

            double G = 1.0 / r;
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

    std::cout << std::setprecision(9) << " # common verts: " << triPair.nCommon << ' '
        << "Rel diff: " << (radNumNum - radNumAnl) / radNumAnl << '\n';
}


void testSingularMFIE(const Mesh::Triangle& obsTri, const Mesh::Triangle& srcTri) {
    assert(config.alpha == 0.0); // only MFIE part should be tested here

    double k2 = config.k * config.k;
    vec3d vobs = obsTri.getVerts()[0];
    vec3d vsrc = srcTri.getVerts()[0];

    pair2i pair = std::minmax(obsTri.getIdx(), srcTri.getIdx());
    Mesh::TriPair triPair(pair);

    // Double numeric integration of k^2/2 * 1/R + 1/R^3
    double radNumNum = 0.0;
    for (const auto& [obs, obsWeight] : obsTri.getQuads()) {
        for (const auto& [src, srcWeight] : srcTri.getQuads()) {
            vec3d R = obs-src;
            double r = R.norm(), r2 = r*r, r3 = r*r2;

            vec3d gradG = R/r3 * (0.5*k2*r2 + 1.0);

            radNumNum -= (obs-vobs).dot(obsTri.getNormal().cross(gradG.cross(src-vsrc)))
                * obsWeight * srcWeight;
        }
    }

    // Numeric-analytic integration of k^2/2 * 1/R + 1/R^3
    double radNumAnl = -obsTri.getSingularMFIE(srcTri, triPair, vobs, vsrc);

    //std::cout << std::setprecision(15)
    //    << "NumNum MFIE: " << radNumNum << ' '
    //    << "NumAnl MFIE: " << radNumAnl << ' '
    //    << "Diff: " << radNumNum - radNumAnl << '\n';

    std::cout << std::setprecision(9) << " # common verts: " << triPair.nCommon << ' '
        << "Rel diff: " << (radNumNum - radNumAnl) / radNumAnl << '\n';
}

/*
void testMassMatrix(const Mesh::Triangle& tri) {
    double k2 = config.k * config.k;
    vec3d v0 = tri.getVerts()[0];
    vec3d v1 = tri.getVerts()[1];

    std::cout << std::setprecision(15);

    // Numerical integration
    auto start = Clock::now();

    double radNum = 0;
    for (const auto& [node, weight] : tri.getQuads()) {
        radNum -= (node-v0).dot(node-v1) * weight;
    }

    Time duration_ms = Clock::now() - start;
    std::cout << "Numeric integration time: " << duration_ms.count() << " ms\n";

    // Analytic integration
    start = Clock::now();

    auto [X0, X1, X2] = tri.getVerts();
    double radAnl =
        -(X0.squaredNorm() + X1.squaredNorm() + X2.squaredNorm() +
            X0.dot(X1) + X0.dot(X2) + X1.dot(X2)) / 12.0
        + (v0 + v1).dot(tri.getCenter()) / 2.0
        - v0.dot(v1) / 2.0;

    duration_ms = Clock::now() - start;
    std::cout << "Analytic integration time: " << duration_ms.count() << " ms\n";

    std::cout
        << "Num: " << radNum << '\n'
        << "Anl: " << radAnl << '\n'
        << "Diff: " << radNum - radAnl << '\n';
}
*/

int main() {
    // ===================== Build sources ==================== //
    std::cout << " Building sources...\n";

    const auto Einc = Exc::importPlaneWaves("config/pwave.txt");
    const auto srcs = importSources(Einc);
    size_t nsrcs = srcs.size();

    // ==================== Build nodes ==================== //
    std::cout << " Building FMM tree...\n";

    auto start0 = Clock::now();
    bool isRootLeaf = nsrcs <= config.maxNodeSrcs;
    auto root = std::make_shared<Node>(srcs, 0, nullptr, isRootLeaf);

    root->buildLists();

    std::cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    std::cout << "   # Leaves: " << leaves.size() << '\n';
    std::cout << "   Max node level: " << Node::getMaxLvl() << "\n\n";

    // ==================== Build nearfield ===================== //
    std::cout << " Building nearfield matrix...     ";

    auto start = Clock::now();
    auto nf = std::make_shared<Nearfield>();
    Time duration_ms = Clock::now() - start;
    std::cout << " in " << duration_ms.count() << " ms\n\n";

    for (const auto& triMap : Mesh::glTriPairs) {
        const auto& [obsTri, srcTri] = triMap.second.getTriPair();

        testSingularMFIE(obsTri, srcTri);
        testSingularMFIE(srcTri, obsTri);
    }

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

