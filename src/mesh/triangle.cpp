#include "triangle.h"

int Mesh::Triangle::numQuads;
std::vector<quadPair<vec3d>> Mesh::Triangle::quadCoeffs;
// std::vector<quadPair<double>> Mesh::Triangle::linQuads;

/* Triangle(iVerts)
 * Construct triangle from global indices of vertices
 * iVerts : global indices of vertices
 */
Mesh::Triangle::Triangle(const vec3i& iVerts, int iTri)
    : iVerts(iVerts), iTri(iTri)
{
    const auto& Xs = getVerts();
    center = (Xs[0] + Xs[1] + Xs[2]) / 3.0;

    for (int i = 0; i < 3; ++i) Ds[i] = Xs[(i+1)%3] - Xs[i];

    nhat = (Ds[0].cross(Ds[1])).normalized();

    area = (Ds[0].cross(-Ds[2])).norm() / 2.0;

    buildTriQuads();
    // buildSelfIntegrated();
    //std::cout << "Built triangle #" << iTri << " with edge lengths " 
    //    << Ds[0].norm() << ", " << Ds[1].norm() << ", " << Ds[2].norm() << '\n';
}

// Return global indices of vertices shared by this and other triangle
int Mesh::Triangle::getNumCommonVerts(const Triangle& other) const {
    if (iTri == other.iTri) return 3;

    int numVerts = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (iVerts[i] == other.iVerts[j])
                ++numVerts;

    assert(numVerts < 3);
    return numVerts;
}

void Mesh::Triangle::buildSelfIntegratedInvR() {
    const auto& [V0, V1, V2] = getVerts();
    double a00 = V0.dot(V0), a01 = V0.dot(V1), a02 = V0.dot(V2);
    double a11 = V1.dot(V1), a12 = V1.dot(V2), a22 = V2.dot(V2);
    double l0 = (V1-V2).norm(), l1 = (V2-V0).norm(), l2 = (V0-V1).norm();

    auto logN = [](double l0, double l1, double l2) {
        double lsum = l0+l1, ldiff = l2-l0;
        return std::log((lsum*lsum - l2*l2) / (l1*l1 - ldiff*ldiff));
        };
    double log0 = logN(l0, l1, l2), log1 = logN(l1, l2, l0), log2 = logN(l2, l0, l1);

    // Integral of \lambda_i * \lambda_i' / r
    auto f0 = [](double l0, double l1, double l2, double log0, double log1, double log2) {
        double l0sq = l0*l0, l1sq = l1*l1, l2sq = l2*l2;
        return log0 / (20.0*l0)
            + log1 * (l0sq + 5.0*l1sq - l2sq) / (120.0*l1sq*l1)
            + log2 * (l0sq - l1sq + 5.0*l2sq) / (120.0*l2sq*l2)
            + (l2-l0) / (60.0*l1sq) + (l1-l0) / (60.0*l2sq);
        };

    // Integral of \lambda_i * \lambda_j' / r
    auto f1 = [](double l0, double l1, double l2, double log0, double log1, double log2) {
        double l0sq = l0*l0, l1sq = l1*l1, l2sq = l2*l2;
        return log2 / (40.0*l2)
            + log0 * (3.0*l0sq + l1sq - l2sq) / (80.0*l0sq*l0)
            + log1 * (l0sq + 3.0*l1sq - l2sq) / (80.0*l1sq*l1)
            + (l2-l1) / (40.0*l0sq) + (l2-l0) / (40.0*l1sq);
        };

    // Integral of \lambda_i / r or \lambda_i' / r
    auto f2 = [](double l0, double l1, double l2, double log0, double log1, double log2) {
        double l0sq = l0*l0, l1sq = l1*l1, l2sq = l2*l2;
        return log0 / (8.0*l0)
            + log1 * (l0sq + 5.0*l1sq - l2sq) / (48.0*l1sq*l1)
            + log2 * (l0sq - l1sq + 5.0*l2sq) / (48.0*l2sq*l2)
            + (l2-l0) / (24.0*l1sq) + (l1-l0) / (24.0*l2sq);
        };

    selfInts[0]
        = f0(l0, l1, l2, log0, log1, log2) * (a00 - 2.0*a01 + a11)
        + f0(l1, l2, l0, log1, log2, log0) * (a00 - 2.0*a02 + a22)
        + f1(l0, l1, l2, log0, log1, log2) * (a00 - a02 - a01 + a12)
        + f1(l0, l1, l2, log0, log1, log2) * (a00 - a01 - a02 + a12);
    selfInts[1] = f2(l0, l1, l2, log0, log1, log2);
    selfInts[2] = f2(l1, l2, l0, log1, log2, log0);
    selfInts[3] = (log0/l0 + log1/l1 + log2/l2) / 3.0;
}

double Mesh::Triangle::getDoubleSelfIntegratedInvR(const vec3d& obsNC, const vec3d& srcNC) const
{
    const auto [V0, V1, V2] = getVerts();
    double a00 = V0.dot(V0), a01 = V0.dot(V1), a02 = V0.dot(V2); // cache?
    const vec3d& sumNC = srcNC + obsNC;
    
    return selfInts[0]
        + selfInts[1] * (-2.0*a00 + 2.0*a01 + (V0-V1).dot(sumNC))
        + selfInts[2] * (-2.0*a00 + 2.0*a02 + (V0-V2).dot(sumNC))
        + selfInts[3] * (a00 - V0.dot(sumNC) + srcNC.dot(obsNC) - 4.0/(k*k));
}

std::pair<double, vec3d>
Mesh::Triangle::getIntegratedInvR(const vec3d& obs, bool doNumeric) const
{
    using namespace Math;

    double scaRad = 0.0;
    vec3d vecRad = vec3d::Zero();
    const vec3d& obsProj = proj(obs);

    if (doNumeric) {
        for (const auto& [node, weight] : triQuads) {
            double dr = (obs-node).norm();
            scaRad += weight / dr;
            vecRad -= weight * (obsProj-proj(node)) / dr;
        }
        return std::make_pair(scaRad, vecRad);
    }

    const auto& Xs = getVerts();
    double d = std::fabs(nhat.dot(obs-Xs[0])), dsq = d*d;
    std::array<vec3d, 3> Ps =
        { proj(Xs[0])-obsProj, proj(Xs[1])-obsProj, proj(Xs[2])-obsProj };

    for (int i = 0; i < 3; ++i) {
        const vec3d& P0 = Ps[i];
        const vec3d& P1 = Ps[(i+1)%3];

        const vec3d& lhat = (P1-P0).normalized();
        const vec3d& uhat = lhat.cross(nhat);
        double l0 = P0.dot(lhat), l1 = P1.dot(lhat);

        const vec3d& P = P0 - l0*lhat;
        double p = P.norm(), p0 = P0.norm(), p1 = P1.norm();
        const vec3d& phat = (fzero(p) ? zeroVec : P.normalized());

        double rsq = p*p + dsq, r0 = std::sqrt(p0*p0 + dsq), r1 = std::sqrt(p1*p1 + dsq);
  
        double f2 = (fzero(r0+l0) || fzero(r1+l1) ?  // use && ?
                std::log(l0/l1) :
                std::log((r1+l1)/(r0+l0)));

        scaRad += phat.dot(uhat) * (
            p*f2 - d * (atan2(p*l1, rsq+d*r1) - atan2(p*l0, rsq+d*r0)));
        vecRad += uhat * (rsq*f2 + l1*r1 - l0*r0);
    }
    scaRad /= 2.0*area;
    vecRad /= 4.0*area;

    return std::make_pair(scaRad, vecRad);
}

std::pair<double, vec3d>
Mesh::Triangle::getIntegratedInvRcubed(const vec3d& obs, bool doNumeric) const
{
    using namespace Math;

    double scaRad3 = 0.0;
    vec3d vecRad3 = vec3d::Zero();
    const vec3d& obsProj = proj(obs);

    if (doNumeric) {
        for (const auto& [node, weight] : triQuads) {
            double dr3 = pow((obs-node).norm(), 3);
            scaRad3 += weight / dr3;
            vecRad3 -= weight * (obsProj-proj(node)) / dr3;
        }
        return std::make_pair(scaRad3, vecRad3);
    }

    const auto& Xs = getVerts();
    double d = std::fabs(nhat.dot(obs-Xs[0])), dsq = d*d;
    std::cout << "d = " << d << '\n';

    std::array<vec3d, 3> Ps =
    { proj(Xs[0])-obsProj, proj(Xs[1])-obsProj, proj(Xs[2])-obsProj };

    for (int i = 0; i < 3; ++i) {
        const vec3d& P0 = Ps[i];
        const vec3d& P1 = Ps[(i+1)%3];

        const vec3d& lhat = (P1-P0).normalized();
        const vec3d& uhat = lhat.cross(nhat);
        double l0 = P0.dot(lhat), l1 = P1.dot(lhat);

        const vec3d& P = P0 - l0*lhat;
        double p = P.norm(), p0 = P0.norm(), p1 = P1.norm();

        const vec3d& phat = (fzero(p) ? zeroVec : P.normalized());

        double rsq = p*p + dsq, r0 = std::sqrt(p0*p0 + dsq), r1 = std::sqrt(p1*p1 + dsq);

        scaRad3 -= phat.dot(uhat) *
            (fzero(d) ?
                (fzero(p) ? 0.0 : l1/(p*r1) - l0/(p*r0)) :
                (atan2(d*l1, p*r1) - atan2(d*l0, p*r0) + atan2(l0, p) - atan2(l1, p)) / d);
        vecRad3 -= uhat *
            (fzero(r0+l0) || fzero(r1+l1) ?  // use && ?
                std::log(l0/l1) :
                std::log((r1+l1)/(r0+l0)));
    }
    scaRad3 /= 2.0*area;
    vecRad3 /= 2.0*area;

    return std::make_pair(scaRad3, vecRad3);
}

double Mesh::Triangle::getDoubleIntegratedSingularEFIE(
    const Triangle& srcTri, const vec3d& obsNC, const vec3d& srcNC) const
{
    const vec3d& srcNCproj = srcTri.proj(srcNC);
    double rad = 0.0;

    for (const auto& [obs, obsWeight] : triQuads) {
        const vec3d& obsProj = srcTri.proj(obs);

        const auto& [scaRad, vecRad] = srcTri.getIntegratedInvR(obs);
        rad += ((obs-obsNC).dot(vecRad+(obsProj-srcNCproj)*scaRad) - 4.0/(k*k)*scaRad)
            * obsWeight;
    }

    return rad;
}

double Mesh::Triangle::getDoubleIntegratedSingularMFIE(
    const Triangle& srcTri, const vec3d& obsNC, const vec3d& srcNC) const
{
    const vec3d& srcNCproj = srcTri.proj(srcNC);
    double rad = 0.0;

    for (const auto& [obs, obsWeight] : triQuads) {
        const vec3d& R = obs-srcNC;
        const vec3d& obsProj = srcTri.proj(obs);

        const auto& [scaRad, vecRad] = srcTri.getIntegratedInvR(obs);
        const auto& [scaRad3, vecRad3] = srcTri.getIntegratedInvRcubed(obs);

        // double check signs and factors
        rad += k*k/2.0
            * (obs-obsNC).dot(R.cross(vecRad+(obsProj-srcNCproj)*scaRad))
            * obsWeight;
        rad +=
            (obs-obsNC).dot(R.cross(vecRad3+(obsProj-srcNCproj)*scaRad3))
            * obsWeight;

        rad += ((obs-obsNC).dot(R.cross(
            vecRad+(obsProj-srcNCproj)*scaRad)))
            * obsWeight;
    }

    return rad;
}

std::pair<cmplx, vec3cd>
Mesh::Triangle::getIntegratedPlaneWave(const vec3d& kvec) const
{
    using namespace Math;

    const auto& Xs = getVerts();
    double alpha = kvec.dot(Ds[0]), beta = -kvec.dot(Ds[2]), gamma = alpha-beta;
    double alphasq = alpha*alpha, betasq = beta*beta;
    cmplx expI_alpha = exp(iu*alpha), expI_beta = exp(iu*beta);
    cmplx // TODO: Only compute if gamma != 0
        f0_alpha = (fzero(alpha) ? -iu : (1.0 - expI_alpha) / alpha),
        f0_beta = (fzero(beta) ? -iu : (1.0 - expI_beta) / beta);
    cmplx
        f1_alpha = (fzero(alpha) ? -0.5 : (1.0 - (1.0 - iu*alpha) * expI_alpha) / alphasq),
        f1_beta = (fzero(beta) ? -0.5 : (1.0 - (1.0 - iu*beta) * expI_beta) / betasq);

    cmplx scaRad; vec3cd vecRad;
    if (fzero(gamma)) {
        cmplx f2 = (fzero(alpha) ? iu/6.0 :
                (expI_alpha*(alphasq + 2.0*iu*alpha - 2.0) + 2.0) / (2.0*alpha*alphasq));
        scaRad = -f1_alpha;
        vecRad = -iu*f2 * (Ds[0] - Ds[2]);
        // radVec = -f1_alpha * (Xs[0] - Xnc[iTri]) - iu*f2 * (Ds[0] - Ds[2]);
    } else {
        cmplx
            I0 = (f0_alpha - f0_beta) / gamma,
            I1 = iu * (I0 + f1_alpha),
            I2 = -iu * (I0 + f1_beta);
        scaRad = I0;
        vecRad = (I1*Ds[0] - I2*Ds[2]) / gamma;
        // radVec = I0 * (Xs[0] - Xnc[iTri]) + (I1*Ds[0] - I2*Ds[2]) / gamma;
    }

    return std::make_pair(scaRad, vecRad);
}

/*
cmplx Mesh::Triangle::getDuffyIntegrated(
    const vec3d& src, const vec3d& obsNC, const vec3d& srcNC) const 
{
    const double k = FMM::wavenum;

    cmplx rad = 0.0;
    // double areaSum = 0.0;
    for (int i = 0; i < 3; ++i) {
        const auto& V0 = glVerts[iVerts[i]];
        const auto& V1 = glVerts[iVerts[(i+1)%3]];
        const double area = (V0-src).cross(V1-src).norm()/2.0; // cache?
        if (Math::fzero(area)) continue;

        // areaSum += area;

        cmplx triRad = 0.0;
        for (const auto& [u, uweight] : linQuads) {
            for (const auto& [v, vweight] : linQuads) {
                const vec3d& obs = src + u*(V0 - src) + u*v*(V1 - V0);

                const double dr = (obs-src).norm();
                const double alpha = dr / u;

                triRad += 
                    ((obs-obsNC).dot(src-srcNC) - 4.0 / (k*k))
                    * exp(iu*k*dr) / alpha // u * Math::helmholtzG(dr, k)
                    * uweight * vweight;
            }
        }

        rad += area * triRad;
    }

    // assert(Math::fzero(areaSum-area));
    // std::cout << areaSum << ' ' << area << '\n';

    return rad;
}*/

/*
void Mesh::Triangle::buildRadMoments() {
    // TODO: Only loop over triangles in nearfield of each other
    for (size_t iObs = 0; iObs < glTris.size(); ++iObs) {
        for (size_t iSrc = 0; iSrc <= iObs; ++iSrc) {
            const auto& obsTri = glTris[iObs];
            const auto& srcTri = glTris[iSrc];

            const auto key = makeUnordered(obsTri.iTri, srcTri.iTri);

            QuadMoments moments = { 0.0, vec3cd::Zero(), vec3cd::Zero(), 0.0 };

            for (const auto& [obs, obsWeight] : obsTri.quads) {
                for (const auto& [src, srcWeight] : srcTri.quads) {
                    const cmplx coeff = 
                        Math::helmholtzG((obs-src).norm(), FMM::wavenum) * obsWeight * srcWeight;
                    get<0>(moments) += obs.dot(src) * coeff;
                    get<1>(moments) += src * coeff;
                    get<2>(moments) += obs * coeff;
                    get<3>(moments) += coeff;
                }
            }

            glRadMoments.emplace(key, moments);
        }
    }
}*/

/*
cmplx Mesh::Triangle::getSurfaceCurrent() const {
    auto triToRWG = triToRWGs[iTri];

    cmplx J = 0.0;
    for (int i = 0; i < 3; ++i) {
        int iRWG = triToRWG.iRWGs[i];
        int isMinus = triToRWG.isMinus[i];
        auto rwg = glSrcs[iRWG];
        const auto& Xnc = rwg->getVertsNC();

        double rwgFunc =
            Math::sign(isMinus) * rwg->leng / (2.0*area) * (center - Xnc[isMinus]);

        J += states.currents[iRWG];
    }
    return J;
}
*/

