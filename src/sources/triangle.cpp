#include "triangle.h"

std::vector<vec3d> Triangle::glVerts;

std::vector<quadPair> Triangle::quadCoeffs;
Precision Triangle::quadPrec;
int Triangle::numQuads;

void Triangle::importVertices(const std::filesystem::path& vpath) {
    std::ifstream file(vpath);
    if (!file) throw std::runtime_error("Unable to find file");
    std::string line;

    while (getline(file,line)) {
        std::istringstream iss(line);
        vec3d vertex;

        if (iss >> vertex)
            glVerts.push_back(vertex);
        else
            throw std::runtime_error("Unable to parse line");
    }
}

TriVec Triangle::importTriangles(
    const std::filesystem::path& vpath, 
    const std::filesystem::path& fpath, 
    Precision prec)
{
    importVertices(vpath);

    quadPrec = prec;
    buildQuadCoeffs(); // build quad coeffs

    std::ifstream file(fpath);
    std::string line;
    if (!file) throw std::runtime_error("Unable to find file");
    TriVec triangles;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        vec3i glIdxs;

        if (iss >> glIdxs)
            triangles.push_back(std::make_shared<Triangle>(glIdxs));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return triangles;
}

/* Triangle(glIdxs)
 * Construct triangle from global indices of vertices
 * glIdxs : global indices of vertices
 */
Triangle::Triangle(const vec3i& glIdxs)
    : glIdxs(glIdxs),
    Xs({ glVerts[glIdxs[0]], glVerts[glIdxs[1]], glVerts[glIdxs[2]] }),
    center( (Xs[0]+Xs[1]+Xs[2])/3.0 )
{
    for (int i = 0; i < 3; ++i) {
        const int ipp = (i+1) % 3;
        Ds[i] = Xs[ipp] - Xs[i];
    }

    nhat = (Ds[0].cross(Ds[1])).normalized();

    // If nhat is pointing inward, reverse it
    // Assume a closed, star-shaped mesh centered at and enclosing the origin
    if (center.dot(nhat) < 0.0) nhat *= -1.0;
    // std::cout << center.dot(nhat) << '\n';

    buildQuads();
};

/* Triangle(idx, X1, X2)
 * Construct triangle on refined mesh from idx of vertex on coarse mesh
 * and midpoint and center of parent triangle
 * idx : index of vertex on coarse mesh
 * X1  : vertex on fine mesh (midpoint or center of parent)
 * X2  : vertex on fine mesh (midpoint or center of parent)
 */
Triangle::Triangle(int idx, const vec3d& X1, const vec3d& X2)
    : Xs( {glVerts[idx], X1, X2} ),
      glIdxs( {idx, 0, 0} ), // TODO: Discard dummy indices
      center((Xs[0]+Xs[1]+Xs[2])/3.0)
{
    for (int i = 0; i < 3; ++i) {
        const int ipp = (i+1) % 3;
        Ds[i] = Xs[ipp] - Xs[i];
    }

    nhat = (Ds[0].cross(Ds[1])).normalized();

    // If nhat is pointing inward, reverse it
    // Assume a closed, star-shaped mesh centered at and enclosing the origin
    if (center.dot(nhat) < 0.0) nhat *= -1.0;

    alpha = acos(
        (Ds[0].normalized()).dot(-Ds[2].normalized())
    );
};

// Construct all 6 sub-tris of this tri
TriArr6 Triangle::getSubtris(const vec3i& glIdxc) {
    TriArr6 subtris;
    int iSubtri = 0;

    for (int i = 0; i < 3; ++i) {
        const int ipp = (i+1) % 3;
        const int idx_i = glIdxc[i], idx_ipp = glIdxc[ipp];
        const auto& mid = (glVerts[idx_i]+glVerts[idx_ipp])/2; // midpoint of ith edge

        for (int j = 0; j < 2; ++j) {
            auto subtri = (!j ?
                std::make_shared<Triangle>(idx_i, mid, center) :
                std::make_shared<Triangle>(idx_ipp, center, mid)
                );

            assert(Math::vecEquals(nhat, subtri->nhat)); // Check orientation
            // std::cout << nhat << '\n' << subtri->nhat << "\n\n";

            subtris[iSubtri] = std::move(subtri);

            ++iSubtri;
        }
    }

    //for (const auto& tri: subtris) {
    //    for (const auto& X : tri->Xs)
    //        std::cout << X << ' ';
    //    std::cout << '\n';
    //}

    return subtris;
}

void Triangle::buildQuadCoeffs() {
    numQuads = [&]() {
        switch (quadPrec) {
            case Precision::VERYLOW: return 1;
            case Precision::LOW:     return 3;
            case Precision::MEDIUM:  return 7;
            case Precision::HIGH:    return 13;
        };
        } ();

    // TODO: Find out why precomputing yields larger error
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
}

//
void Triangle::buildQuads() {
    auto baryToPos = [&](const vec3d& ws) {
        return ws[0]*Xs[0] + ws[1]*Xs[1] + ws[2]*Xs[2];
    };

    quads.reserve(numQuads);

    for (const auto& [coeffs, weight] : quadCoeffs) {
        quads.emplace_back(baryToPos(coeffs), weight);
    }
}
//

/*
void Triangle::buildQuads() {
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