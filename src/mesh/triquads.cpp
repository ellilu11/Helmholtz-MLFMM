#include "triangle.h"

/* buildQuadCoeffs(prec)
 * Build quadrature nodes and weights for given precision
 * prec : precision level
 */
void Mesh::Triangle::buildQuadCoeffs(Precision prec) {
    numQuads = getNumQuads(prec);
    quadCoeffs.reserve(numQuads);

    switch (prec) {
        // 1-point
        case Precision::VERYLOW: { 
            constexpr double weight = 1.0/2.0;
            const vec3d ws0(1.0/3.0, 1.0/3.0, 1.0/3.0);
            quadCoeffs.emplace_back(ws0, weight);
            break;
        }

        // 3-point
        case Precision::LOW: { 
            constexpr double weight0 = 1.0/6.0;
            std::vector<vec3d> wss0;
            wss0.emplace_back(2.0/3.0, 1.0/6.0, 1.0/6.0);
            wss0.emplace_back(1.0/6.0, 2.0/3.0, 1.0/6.0);
            wss0.emplace_back(1.0/6.0, 1.0/6.0, 2.0/3.0);
            for (const auto& ws : wss0)
                quadCoeffs.emplace_back(ws, weight0);
            break;
        }
    
        // 4-point
        case Precision::MEDLOW: { 
            constexpr double weight0 = -0.281250000000000;
            const vec3d ws0(1.0/3.0, 1.0/3.0, 1.0/3.0);
            quadCoeffs.emplace_back(ws0, weight0);

            constexpr double weight1 = 0.260416666666667;
            constexpr double alpha1 = 0.600000000000000, beta1 = 0.200000000000000;
            std::vector<vec3d> wss1;
            wss1.emplace_back(alpha1, beta1, beta1);
            wss1.emplace_back(beta1, alpha1, beta1);
            wss1.emplace_back(beta1, beta1, alpha1);
            for (const auto& ws : wss1)
                quadCoeffs.emplace_back(ws, weight1);

            constexpr double weightErr = weight0 + 3.0*(weight1) - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }

        // 6-point
        case Precision::MEDIUM: {
            constexpr double weight1 = 0.111690794839005;
            constexpr double alpha1 = 0.108103018168070, beta1 = 0.445948490915965;
            std::vector<vec3d> wss1;
            wss1.emplace_back(alpha1, beta1, beta1);
            wss1.emplace_back(beta1, alpha1, beta1);
            wss1.emplace_back(beta1, beta1, alpha1);
            for (const auto& ws : wss1)
                quadCoeffs.emplace_back(ws, weight1);

            constexpr double weight2 = 0.054975871827661;
            constexpr double alpha2 = 0.816847572980459, beta2 = 0.091576213509771;
            std::vector<vec3d> wss2;
            wss2.emplace_back(alpha2, beta2, beta2);
            wss2.emplace_back(beta2, alpha2, beta2);
            wss2.emplace_back(beta2, beta2, alpha2);
            for (const auto& ws : wss2)
                quadCoeffs.emplace_back(ws, weight2);

            constexpr double weightErr = 3.0*(weight1+weight2) - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }
    
        // 7-point
        case Precision::MEDHIGH: { 
            constexpr double weight0 = 0.1125;
            const vec3d ws0(1.0/3.0, 1.0/3.0, 1.0/3.0);
            quadCoeffs.emplace_back(ws0, weight0);

            constexpr double weight1 = 0.066197076394253;
            constexpr double alpha1 = 0.059715871789770, beta1 = 0.470142064105115;
            std::vector<vec3d> wss1;
            wss1.emplace_back(alpha1, beta1, beta1);
            wss1.emplace_back(beta1, alpha1, beta1);
            wss1.emplace_back(beta1, beta1, alpha1);
            for (const auto& ws : wss1)
                quadCoeffs.emplace_back(ws, weight1);

            constexpr double weight2 = 0.0629695902724135;
            constexpr double alpha2 = 0.797426985353087, beta2 = 0.101286507323456;
            std::vector<vec3d> wss2;
            wss2.emplace_back(alpha2, beta2, beta2);
            wss2.emplace_back(beta2, alpha2, beta2);
            wss2.emplace_back(beta2, beta2, alpha2);
            for (const auto& ws : wss2)
                quadCoeffs.emplace_back(ws, weight2);

            constexpr double weightErr = weight0 + 3.0*(weight1+weight2) - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }

        // 12-point
        case Precision::HIGH: { 
            constexpr double weight1 = 0.058393137863189;
            constexpr double alpha1 = 0.501426509658179, beta1 = 0.249286745170910;
            vec3d ws1(alpha1, beta1, beta1);
            std::vector<vec3d> wss1;
            wss1.emplace_back(alpha1, beta1, beta1);
            wss1.emplace_back(beta1, alpha1, beta1);
            wss1.emplace_back(beta1, beta1, alpha1);
            for (const auto& ws : wss1)
                quadCoeffs.emplace_back(ws, weight1);

            constexpr double weight2 = 0.025422453185103;
            constexpr double alpha2 = 0.873821971016996, beta2 = 0.063089014491502;
            vec3d ws2(alpha2, beta2, beta2);
            std::vector<vec3d> wss2;
            wss2.emplace_back(alpha2, beta2, beta2);
            wss2.emplace_back(beta2, alpha2, beta2);
            wss2.emplace_back(beta2, beta2, alpha2);
            for (const auto& ws : wss2)
                quadCoeffs.emplace_back(ws, weight2);

            constexpr double weight3 = 0.041425537809187;
            constexpr double alpha3 = 0.053145049844817, beta3 = 0.310352451033784, gamma = 0.636502499121399;
            vec3d ws3(alpha3, beta3, gamma);
            std::vector<vec3d> wss3;
            wss3.emplace_back(alpha3, beta3, gamma);
            wss3.emplace_back(beta3, gamma, alpha3);
            wss3.emplace_back(gamma, alpha3, beta3);

            wss3.emplace_back(beta3, alpha3, gamma);
            wss3.emplace_back(alpha3, gamma, beta3);
            wss3.emplace_back(gamma, beta3, alpha3);
            for (const auto& ws : wss3)
                quadCoeffs.emplace_back(ws, weight3);

            constexpr double weightErr = 3.0*(weight1+weight2) + 6.0*weight3 - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }

        // 13-point
        case Precision::VERYHIGH: { 
            constexpr double weight0 = -0.074785022233841;
            const vec3d ws0(1.0/3.0, 1.0/3.0, 1.0/3.0);
            quadCoeffs.emplace_back(ws0, weight0);

            constexpr double weight1 = 0.087807628716604;
            constexpr double alpha1 = 0.479308067841920, beta1 = 0.260345966079040;
            vec3d ws1(alpha1, beta1, beta1);
            std::vector<vec3d> wss1;
            wss1.emplace_back(alpha1, beta1, beta1);
            wss1.emplace_back(beta1, alpha1, beta1);
            wss1.emplace_back(beta1, beta1, alpha1);
            for (const auto& ws : wss1)
                quadCoeffs.emplace_back(ws, weight1);

            constexpr double weight2 = 0.026673617804419;
            constexpr double alpha2 = 0.869739794195568, beta2 = 0.065130102902216;
            vec3d ws2(alpha2, beta2, beta2);
            std::vector<vec3d> wss2;
            wss2.emplace_back(alpha2, beta2, beta2);
            wss2.emplace_back(beta2, alpha2, beta2);
            wss2.emplace_back(beta2, beta2, alpha2);
            for (const auto& ws : wss2)
                quadCoeffs.emplace_back(ws, weight2);

            constexpr double weight3 = 0.038556880445128;
            constexpr double alpha3 = 0.048690315425316, beta3 = 0.312865496004874, gamma = 0.638444188569810;
            vec3d ws3(alpha3, beta3, gamma);
            std::vector<vec3d> wss3;
            wss3.emplace_back(alpha3, beta3, gamma);
            wss3.emplace_back(beta3, gamma, alpha3);
            wss3.emplace_back(gamma, alpha3, beta3);

            wss3.emplace_back(beta3, alpha3, gamma);
            wss3.emplace_back(alpha3, gamma, beta3);
            wss3.emplace_back(gamma, beta3, alpha3);
            for (const auto& ws : wss3)
                quadCoeffs.emplace_back(ws, weight3);

            constexpr double weightErr = weight0 + 3.0*(weight1+weight2) + 6.0*weight3 - 0.5;
            static_assert(weightErr > -Math::FEPS && weightErr < Math::FEPS);
            break;
        }

        /* // 19-point, for debugging
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
        */
    }

    assert(quadCoeffs.size() == numQuads);
}

