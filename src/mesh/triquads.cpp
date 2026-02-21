#include "triangle.h"

void Mesh::Triangle::buildQuadCoeffs(Precision prec) {
    numQuads = [&]() {
        switch (prec) {
            case Precision::VERYLOW:  return 1;
            case Precision::LOW:      return 3;
            case Precision::MEDIUM:   return 7;
            case Precision::HIGH:     return 13;
            case Precision::VERYHIGH: return 19;
        };
        } ();

    // TODO: Find out why precomputing yields larger error
    quadCoeffs.reserve(numQuads);

    switch (prec) {
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

        case Precision::VERYHIGH: {
            constexpr double weight0 = 0.048567898141400;
            const vec3d ws0(1.0/3.0, 1.0/3.0, 1.0/3.0);
            quadCoeffs.emplace_back(ws0, weight0);

            constexpr double weight1 = 0.015667350113570;
            constexpr double alpha1 = 0.020634961602525, beta1 = 0.489682519198738;
            vec3d ws1(alpha1, beta1, beta1);
            std::vector<vec3d> wss1;
            Math::buildPermutations(ws1, wss1, 0);
            for (const auto& ws : wss1)
                quadCoeffs.emplace_back(ws, weight1);

            constexpr double weight2 = 0.038913770502387;
            constexpr double alpha2 = 0.125820817014127, beta2 = 0.437089591492937;
            vec3d ws2(alpha2, beta2, beta2);
            std::vector<vec3d> wss2;
            Math::buildPermutations(ws2, wss2, 0);
            for (const auto& ws : wss2)
                quadCoeffs.emplace_back(ws, weight2);

            constexpr double weight3 = 0.039823869463605;
            constexpr double alpha3 = 0.623592928761935, beta3 = 0.188203535619033;
            vec3d ws3(alpha3, beta3, beta3);
            std::vector<vec3d> wss3;
            Math::buildPermutations(ws3, wss3, 0);
            for (const auto& ws : wss3)
                quadCoeffs.emplace_back(ws, weight3);

            constexpr double weight4 = 0.012788837829349;
            constexpr double alpha4 = 0.910540973211095, beta4 = 0.044729513394453;
            vec3d ws4(alpha4, beta4, beta4);
            std::vector<vec3d> wss4;
            Math::buildPermutations(ws4, wss4, 0);
            for (const auto& ws : wss4)
                quadCoeffs.emplace_back(ws, weight4);

            constexpr double weight5 = 0.021641769688645;
            constexpr double alpha5 = 0.036838412054736, beta5 = 0.221962989160766, gamma = 0.741198598784498;
            vec3d ws5(alpha5, beta5, gamma);
            std::vector<vec3d> wss5;
            Math::buildPermutations(ws5, wss5, 0);
            for (const auto& ws : wss5)
                quadCoeffs.emplace_back(ws, weight5);

            constexpr double weightErr = weight0 + 3.0*(weight1+weight2+weight3+weight4) + 6.0*weight5 - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }
    }

    assert(quadCoeffs.size() == numQuads);

    /* Also build linear quadrature nodes and weights (for Duffy integration)
    constexpr int numLinQuads = 3;
    auto [nodes, weights] = Math::gaussLegendre(numLinQuads, 0.0, 1.0);
    assert(nodes.size() == numLinQuads);
    for (int i = 0; i < numLinQuads; ++i)
        linQuads.emplace_back(nodes[i], weights[i]);
    */
}

// Build triangular quadrature nodes and weights
void Mesh::Triangle::buildTriQuads() {
    const auto& Xs = getVerts();

    auto baryToPos = [&](const vec3d& ws) {
        return ws[0]*Xs[0] + ws[1]*Xs[1] + ws[2]*Xs[2];
    };

    triQuads.reserve(numQuads);
    for (const auto& [coeffs, weight] : quadCoeffs)
        triQuads.emplace_back(baryToPos(coeffs), weight);
}

/*
void Mesh::Triangle::buildQuads() {
    auto baryToPos = [&](double w0, double w1, double w2) {
        return w0*Xs[0] + w1*Xs[1] + w2*Xs[2];
    };

    quads.reserve(numQuads);

    switch (quadPrec) {
        case Precision::VERYLOW:
            quads.emplace_back(center, 1.0/2.0);
            break;

        case Precision::LOW: {
            constexpr double weight = 1.0/6.0;
            quads.emplace_back(baryToPos(2.0/3.0, 1.0/6.0, 1.0/6.0), weight);
            quads.emplace_back(baryToPos(1.0/6.0, 2.0/3.0, 1.0/6.0), weight);
            quads.emplace_back(baryToPos(1.0/6.0, 1.0/6.0, 2.0/3.0), weight);
            break;
        }

        case Precision::MEDIUM: {
            constexpr double weight0 = 0.1125;
            quads.emplace_back(center, weight0);

            constexpr double weight1 = 0.066197076394253;
            constexpr double alpha = 0.059715871789770, beta = 0.470142064105115;
            quads.emplace_back(baryToPos(alpha, beta, beta), weight1);
            quads.emplace_back(baryToPos(beta, alpha, beta), weight1);
            quads.emplace_back(baryToPos(beta, beta, alpha), weight1);

            constexpr double weight2 = 0.0629695902724135;
            constexpr double gamma = 0.797426985353087, delta = 0.101286507323456;
            quads.emplace_back(baryToPos(gamma, delta, delta), weight2);
            quads.emplace_back(baryToPos(delta, gamma, delta), weight2);
            quads.emplace_back(baryToPos(delta, delta, gamma), weight2);

            constexpr double weightErr = weight0 + 3.0*(weight1+weight2) - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }

        case Precision::HIGH: {
            constexpr double weight0 = -0.074785022233841;
            quads.emplace_back(center, weight0);

            constexpr double weight1 = 0.087807628716604;
            constexpr double alpha = 0.479308067841920, beta = 0.260345966079040;
            quads.emplace_back(baryToPos(alpha, beta, beta), weight1);
            quads.emplace_back(baryToPos(beta, alpha, beta), weight1);
            quads.emplace_back(baryToPos(beta, beta, alpha), weight1);

            constexpr double weight2 = 0.026673617804419;
            constexpr double gamma = 0.869739794195568, delta = 0.065130102902216;
            quads.emplace_back(baryToPos(gamma, delta, delta), weight2);
            quads.emplace_back(baryToPos(delta, gamma, delta), weight2);
            quads.emplace_back(baryToPos(delta, delta, gamma), weight2);

            constexpr double weight3 = 0.0385568804451285;
            constexpr double eps = 0.048690315425316, zeta = 0.312865496004874, theta = 0.638444188569810;
            quads.emplace_back(baryToPos(eps, zeta, theta), weight3);
            quads.emplace_back(baryToPos(theta, eps, zeta), weight3);
            quads.emplace_back(baryToPos(zeta, theta, eps), weight3);
            quads.emplace_back(baryToPos(eps, theta, zeta), weight3);
            quads.emplace_back(baryToPos(zeta, eps, theta), weight3);
            quads.emplace_back(baryToPos(theta, zeta, eps), weight3);

            constexpr double weightErr = weight0 + 3.0*(weight1+weight2) + 6.0*weight3 - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }
    }

    assert(quads.size() == numQuads);
}
*/