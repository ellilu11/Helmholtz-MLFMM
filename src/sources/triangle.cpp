#include "triangle.h"

std::vector<quadPair> Triangle::quadCoeffs;
int Triangle::numQuads;

Triangle::Triangle(
    const vec3i& vIdx,
    const std::vector<vec3d>& vList,
    const Precision quadPrec)
    : vIdx(vIdx),
    Xs({ vList[vIdx[0]], vList[vIdx[1]], vList[vIdx[2]] }),
    center( (Xs[0]+Xs[1]+Xs[2])/3.0 )
{
    for (int i = 0; i < 3; ++i) {
        const int ipp = Math::wrapIdxToRange(i+1, 3);
        Ds[i] = Xs[ipp] - Xs[i];
    }

    nhat = (Ds[0].cross(Ds[1])).normalized();

    buildQuads(quadPrec);
};

void Triangle::buildQuadCoeffs(Precision quadPrec) {
    numQuads = [&]() {
        switch (quadPrec) {
            case Precision::VERYLOW: return 1;
            case Precision::LOW:     return 3;
            case Precision::MEDIUM:  return 7;
            case Precision::HIGH:    return 13;
        };
        } ();

    quadCoeffs.reserve(numQuads);

    switch (quadPrec) {
        case Precision::VERYLOW: {
            constexpr double weight = 1.0/2.0;
            const vec3d ws0(1.0/3.0, 1.0/3.0, 1.0/3.0);
            quadCoeffs.emplace_back(ws0, weight);
            break;
        }

        case Precision::LOW: {
            constexpr double weight0 = 1.0/6.0;
            // vec3d ws0(2.0/3.0, 1.0/6.0, 1.0/6.0);
            std::vector<vec3d> wss0;
            wss0.emplace_back(2.0/3.0, 1.0/6.0, 1.0/6.0);
            wss0.emplace_back(1.0/6.0, 2.0/3.0, 1.0/6.0);
            wss0.emplace_back(1.0/6.0, 1.0/6.0, 2.0/3.0);
            // Math::buildPermutations(ws0, wss0, 0);
            for (const auto& ws : wss0)
                quadCoeffs.emplace_back(ws, weight0);
            break;
        }

        case Precision::MEDIUM: {
            constexpr double weight0 = 0.1125;
            const vec3d ws0(1.0/3.0, 1.0/3.0, 1.0/3.0);
            quadCoeffs.emplace_back(ws0, weight0);

            constexpr double weight1 = 0.066197076394253;
            constexpr double alpha1 = 0.059715871789770, beta1 = 0.470142064105115;
            // vec3d ws1(alpha1, beta1, beta1);
            std::vector<vec3d> wss1;
            wss1.emplace_back(alpha1, beta1, beta1);
            wss1.emplace_back(beta1, alpha1, beta1);
            wss1.emplace_back(beta1, beta1, alpha1);
            // Math::buildPermutations(ws1, wss1, 0);
            for (const auto& ws : wss1)
                quadCoeffs.emplace_back(ws, weight1);

            constexpr double weight2 = 0.0629695902724135;
            constexpr double alpha2 = 0.797426985353087, beta2 = 0.101286507323456;
            // vec3d ws2(alpha2, beta2, beta2);
            std::vector<vec3d> wss2;
            wss2.emplace_back(alpha2, beta2, beta2);
            wss2.emplace_back(beta2, alpha2, beta2);
            wss2.emplace_back(beta2, beta2, alpha2);
            // Math::buildPermutations(ws2, wss2, 0);
            for (const auto& ws : wss2)
                quadCoeffs.emplace_back(ws, weight2);

            constexpr double weightErr = weight0 + 3.0*(weight1+weight2) - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }

        case Precision::HIGH: { // TODO: Fix buildPermutations before selecting this option
            constexpr double weight0 = -0.074785022233841;
            const vec3d ws0(1.0/3.0, 1.0/3.0, 1.0/3.0);
            quadCoeffs.emplace_back(ws0, weight0);

            constexpr double weight1 = 0.087807628716604;
            constexpr double alpha1 = 0.479308067841920, beta1 = 0.260345966079040;
            vec3d ws1(alpha1, beta1, beta1);
            std::vector<vec3d> wss1;
            Math::buildPermutations(ws1, wss1, 0);
            for (const auto& ws : wss1)
                quadCoeffs.emplace_back(ws, weight1);

            constexpr double weight2 = 0.026673617804419;
            constexpr double alpha2 = 0.869739794195568, beta2 = 0.065130102902216;
            vec3d ws2(alpha2, beta2, beta2);
            std::vector<vec3d> wss2;
            Math::buildPermutations(ws2, wss2, 0);
            for (const auto& ws : wss2)
                quadCoeffs.emplace_back(ws, weight2);

            constexpr double weight3 = 0.0385568804451285;
            constexpr double alpha3 = 0.048690315425316, beta3 = 0.312865496004874, gamma = 0.638444188569810;
            vec3d ws3(alpha3, beta3, gamma);
            std::vector<vec3d> wss3;
            Math::buildPermutations(ws3, wss3, 0);
            for (const auto& ws : wss3)
                quadCoeffs.emplace_back(ws, weight3);

            constexpr double weightErr = weight0 + 3.0*(weight1+weight2) + 6.0*weight3 - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }
    }

    assert(quadCoeffs.size() == numQuads);

    //for (const auto& [coeffs, weight] : quadCoeffs)
    //    std::cout << coeffs << '\n';
}

void Triangle::buildQuads(Precision quadPrec) {
    auto coeffsToNode = [&](const vec3d& ws) {
        return ws[0]*Xs[0] + ws[1]*Xs[1] + ws[2]*Xs[2];
    };
    //auto coeffsToNode = [&](double w0, double w1, double w2) {
    //    return w0*Xs[0] + w1*Xs[1] + w2*Xs[2];
    //};

    quads.reserve(numQuads);

    for (const auto& [coeffs, weight] : quadCoeffs) {
        quads.emplace_back(coeffsToNode(coeffs), weight);
        // quads.emplace_back(coeffsToNode(coeffs[0], coeffs[1], coeffs[2]), weight);
    }

    //for (const auto& [node, weight] : quads)
    //    std::cout << std::setprecision(15) << node << '\n';
    //std::cout << "\n\n";
}