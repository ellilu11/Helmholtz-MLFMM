#include "triangle.h"

std::vector<vec3d> Triangle::glVerts;
std::vector<Triangle> Triangle::glTris;
PairHashMap<int> Triangle::glEdgeToIdx;
PairHashMap<vec2i> Triangle::glEdgeToTri;

size_t Triangle::numCoarseVerts;
size_t Triangle::numCoarseTris;

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

        if (iss >> vertex) glVerts.push_back(vertex);
        else throw std::runtime_error("Unable to parse line");
    }
}

void Triangle::importTriangles(
    const std::filesystem::path& vpath, 
    const std::filesystem::path& fpath, 
    Precision prec)
{
    importVertices(vpath);
    numCoarseVerts = glVerts.size();

    quadPrec = prec;
    buildQuadCoeffs();

    std::ifstream file(fpath);
    std::string line;
    if (!file) throw std::runtime_error("Unable to find file");

    int iTri = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        vec3i iVerts;

        if (iss >> iVerts) glTris.emplace_back(iVerts, iTri++);
        else throw std::runtime_error("Unable to parse line");
    }
    numCoarseTris = glTris.size();

    // TODO: Merge these 3 functions
    refineVertices();
    buildSubtris();
    buildEdgeToTri();
}

/* Triangle(iVerts)
 * Construct triangle from global indices of vertices
 * iVerts : global indices of vertices
 */
Triangle::Triangle(const vec3i& iVerts, int iTri) 
    : iVerts(iVerts), iTri(iTri)
{
    assert(iVerts.maxCoeff() < glVerts.size());
    //std::cout << "Building triangle with indices: " 
    //          << iVerts[0] << ", " 
    //          << iVerts[1] << ", " 
    //          << iVerts[2] << '\n';

    const auto Xs = getVerts();

    for (int i = 0; i < 3; ++i) Ds[i] = Xs[(i+1)%3] - Xs[i];
    nhat = (Ds[0].cross(Ds[1])).normalized();

    // If nhat is pointing inward, reverse it
    // Assume a closed, star-shaped mesh centered at and enclosing the origin
    if (center.dot(nhat) < 0.0) nhat *= -1.0;
    // std::cout << center.dot(nhat) << '\n';

    //alpha = acos(
    //    (Ds[0].normalized()).dot(-Ds[2].normalized())
    //);

    buildQuads(Xs);
}

void Triangle::refineVertices() {
    assert(!glVerts.empty());

    int numVerts = glVerts.size();
    for (auto& tri : glTris) {
        // Add center vertex
        tri.iCenter = numVerts++;
        glVerts.push_back(tri.center); 

        // Add midpoints of edges
        for (int i = 0; i < 3; ++i) {
            const int idx0 = tri.iVerts[i], idx1 = tri.iVerts[(i+1)%3];
            const vec3d mid = (glVerts[idx0] + glVerts[idx1]) / 2.0;
            const auto& edge = makeUnordered(idx0,idx1);

            // Check if midpoint already exists
            if (glEdgeToIdx.find(edge) != glEdgeToIdx.end())
                continue;

            glEdgeToIdx.emplace(edge, numVerts++);
            glVerts.push_back(mid);
        }
    }

    assert(numVerts == glVerts.size());
    // std::cout << "Refined to " << glVerts.size() << " vertices\n";
}

// Construct all 6 subtris of all coarse tris and add to global list
void Triangle::buildSubtris(){
    int iTri = glTris.size();

    std::vector<Triangle> glSubtris;
    glSubtris.reserve(6 * glTris.size());

    for (const auto& tri : glTris) {
        for (int i = 0; i < 3; ++i) {
            const int idx0 = tri.iVerts[i], idx1 = tri.iVerts[(i+1)%3], idxCenter = tri.iCenter;
            const int idxMid = glEdgeToIdx[makeUnordered(idx0, idx1)];

            for (int j = 0; j < 2; ++j) {
                // Assign vertex on coarse mesh as first vertex of subtri
                const auto& iVerts = (!j ?
                    vec3i(idx0, idxMid, idxCenter) :
                    vec3i(idx1, idxCenter, idxMid)
                    );

                glSubtris.emplace_back(iVerts, iTri);
                // assert(Math::vecEquals(nhat,glSubtris[iTri]->nhat)); // Check orientation
                // std::cout << nhat << '\n' << glSubtris[iTri]->nhat << "\n\n";
                ++iTri;
            }
        }
    }

    // Append list of subtris to global list of tris
    glTris.insert(glTris.end(), glSubtris.begin(), glSubtris.end());
}

void Triangle::buildEdgeToTri() {
    const auto glSubtris = glTris | std::views::drop(numCoarseTris);

    for (const auto& tri : glSubtris) {
        for (int i = 0; i < 3; ++i) {
            const int idx0 = tri.iVerts[i], idx1 = tri.iVerts[(i+1)%3];
            const auto& edge = makeUnordered(idx0, idx1);

            if (glEdgeToTri.find(edge) == glEdgeToTri.end())
                glEdgeToTri.emplace(edge, vec2i(tri.iTri, -1));
            else
                glEdgeToTri[edge][1] = tri.iTri;
        }
    }
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
void Triangle::buildQuads(const std::array<vec3d,3>& Xs) {
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