#include "rwg.h"

RWGVec RWG::glSubrwgs;
std::vector<std::vector<int>> RWG::vertsToSubrwgs;

SrcVec RWG::importRWG(
    const std::filesystem::path& vpath,
    const std::filesystem::path& tpath,
    const std::filesystem::path& rpath,
    const Precision quadPrec,
    const std::shared_ptr<Excitation::PlaneWave> Einc)
{
    auto triangles = Triangle::importTriangles(vpath,tpath,quadPrec);

    std::ifstream file(rpath);
    std::string line;
    if (!file) throw std::runtime_error("Unable to find file");
    SrcVec rwgs;
    size_t rwgIdx = 0;

    while (getline(file,line)) {
        std::istringstream iss(line);
        Eigen::Vector4i idx4;

        if (iss >> idx4)
            rwgs.emplace_back(make_shared<RWG>(Einc,rwgIdx++,idx4,triangles));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return rwgs;
}

RWG::RWG(
    std::shared_ptr<Excitation::PlaneWave> Einc,
    size_t rwgIdx,
    const Eigen::Vector4i& idx4,
    const TriVec& triangles)
    : Source(std::move(Einc), rwgIdx),
      tris({triangles[idx4[2]], triangles[idx4[3]]}),
      idx_c({idx4[0],idx4[1]}),
      bc(std::make_unique<BC>(this))
{
    Xc[0] = Triangle::glVerts[idx_c[0]], Xc[1] = Triangle::glVerts[idx_c[1]];
    center = (Xc[0]+Xc[1])/2.0;
    leng = (Xc[0]-Xc[1]).norm();

    // Find non-common vertices
    for (int i = 0; i < 2; ++i)
        for (const auto& glIdxs : tris[i]->glIdxs)
            if (glIdxs != idx_c[0] && glIdxs != idx_c[1]) {
                idx_nc[i] = glIdxs;
                Xnc[i] = Triangle::glVerts[glIdxs];
            }

    buildSubRWGs();

    buildVoltage(); // needs Xnc initialized!

    //std::cout << idx_c[0] << ' '<< idx_c[1] << ' '
    //    << idx_c_pm[0] << ' ' << idx_c_pm[1] << "\n";

    //std::cout << '(' << X[0] << ") (" << X[1] << ") ("
    //    << Xnc[0] << ") (" << Xnc[1] << ") " << leng << '\n';
};

RWG::RWG(
    std::shared_ptr<Triangle> tri0,
    std::shared_ptr<Triangle> tri1)
    : tris( {std::move(tri0), std::move(tri1)} )
{ 
    using namespace Math;

    // Find common vertices
    int k = 0;
    for (int i = 0; i < 3; ++i) {
        const vec3d& X0 = tris[0]->Xs[i];
        for (int j = 0; j < 3; ++j) {
            const vec3d& X1 = tris[1]->Xs[j];
            if (vecEquals(X0,X1)) Xc[k++] = X0;
        }
    }
    const vec3d& dX = Xc[1]-Xc[0];
    leng = dX.norm();

    // Find non-common vertices
    for (int i = 0; i < 2; ++i)
        for (const auto& X : tris[i]->Xs)
            if (!vecEquals(X,Xc[0]) && !vecEquals(X,Xc[1]))
                Xnc[i] = X;

    // Reorder tris as needed
    const vec3d& nhat0 = dX.cross(Xnc[0] - Xc[0]);
    const vec3d& nhat1 = dX.cross(Xnc[1] - Xc[0]);
    assert(nhat0.dot(nhat0 - nhat1) > 0);
    // std::cout << nhat0.dot(nhat0 - nhat1) << '\n';
    // if (nhat0.dot(nhat0 - nhat1) < 0) std::swap(tris[0],tris[1]);
}

void RWG::buildSubRWGs() {
    const Eigen::MatrixXi subIdxs{
        {0, 0, 4, 1, 2, 3},
        {5, 1, 5, 2, 3, 4}
    };

    int iTri = 0;
    TriVec midSubtris;

    for (const auto& tri : tris) {
        const auto& subtris = tri->getSubtris(vec3i(idx_nc[iTri], idx_c[0], idx_c[1]));
        const int iTriBy8 = 8*iTri;
        
        for (int iSub = 0; iSub < 6; ++iSub) {
            auto subrwg = 
                std::make_shared<RWG>(subtris[subIdxs(0,iSub)], subtris[subIdxs(1,iSub)]);

            subrwgs[iTriBy8+iSub] = subrwg; // TODO: Move assign
            glSubrwgs.push_back(std::move(subrwg));
                
            if (iSub == 2 || iSub == 3)
                midSubtris.push_back(subtris[iSub]);
        }

        ++iTri;
    }
    
    // Construct RWGs along common edge of parent RWG
    for (int iSub = 0; iSub < 2; ++iSub) {
        auto subrwg = std::make_shared<RWG>(midSubtris[iSub],midSubtris[iSub+2]);
        subrwgs[iSub+6] = subrwg; // TODO: Move assign
        glSubrwgs.push_back(std::move(subrwg));
    }
    //subrwgs[6] = std::make_shared<RWG>(midSubtris[0], midSubtris[2]);
    //subrwgs[7] = std::make_shared<RWG>(midSubtris[1], midSubtris[3]);

    //for (const auto& rwg : subrwgs)
    //    std::cout << rwg->Xc[0] << ' ' << rwg->Xc[1] << '\n';
}

void RWG::buildVertsToSubrwgs(int numVerts) {
    vertsToSubrwgs.resize(numVerts);

    int iTri = 0;
    for (const auto& rwg : glSubrwgs) {
        // int glIdxs = tri->glIdxs[0]; // only examine vert of subtri in coarse mesh
        // vertsToSubrwgs[glIdxs].push_back(iTri);
        ++iTri;
    }

    //int glIdxs = 0;
    //for (const auto& verts : vertsToSubtris)
    //    std::cout << glIdxs++ << " " << verts.size() << '\n';
}


vec3cd RWG::getIntegratedPlaneWave(const vec3d& kvec, bool doNumeric) const {
    using namespace Math;

    vec3cd rad = vec3cd::Zero();
    int triIdx = 0;

    if (doNumeric) {
        for (const auto& tri : tris) {
            for (const auto& [node, weight] : tri->getQuads())
                rad += weight * exp(iu*kvec.dot(node))
                        * (node - Xnc[triIdx])
                        * sign(triIdx);
            ++triIdx;
        }

        return leng * rad;
    }

    for (const auto& tri : tris) {
        const auto& Xs = tri->Xs, Ds = tri->Ds;

        const double alpha = kvec.dot(Ds[0]), beta = -kvec.dot(Ds[2]), gamma = alpha-beta;

        const double alphasq = alpha*alpha, betasq = beta*beta;

        const cmplx expI_alpha = exp(iu*alpha), expI_beta = exp(iu*beta);

        const cmplx // TODO: Only compute if gamma != 0
            f0_alpha = (approxZero(alpha) ? -iu : (1.0 - expI_alpha) / alpha),
            f0_beta = (approxZero(beta) ? -iu : (1.0 - expI_beta) / beta);

        const cmplx
            f1_alpha = (approxZero(alpha) ? -0.5 : (1.0 - (1.0 - iu*alpha) * expI_alpha) / alphasq),
            f1_beta = (approxZero(beta) ? -0.5 : (1.0 - (1.0 - iu*beta) * expI_beta) / betasq);

        vec3cd radVec;

        if (approxZero(gamma)) {
            const cmplx
                f2 = (expI_alpha*(alphasq + 2.0*iu*alpha - 2.0) + 2.0) / (2.0*alpha*alphasq);

            radVec = -f1_alpha * (Xs[0] - Xnc[triIdx]) - iu*f2 * (Ds[0] - Ds[2]);

        } else {
            const cmplx
                I0 = (f0_alpha - f0_beta) / gamma,
                I1 = iu * (I0 + f1_alpha),
                I2 = -iu * (I0 + f1_beta);

            radVec = I0 * (Xs[0] - Xnc[triIdx]) + (I1*Ds[0] - I2*Ds[2]) / gamma;
        }

        rad += exp(iu*kvec.dot(Xs[0])) * radVec * sign(triIdx++);
    }

    return leng * rad;
}

/* getIntegratedRad(src)
 * Return the radiated field due to src tested with this RWG
 */
cmplx RWG::getIntegratedRad(const std::shared_ptr<Source> src) const {

    const auto srcRWG = dynamic_pointer_cast<RWG>(src);
    const double k = Einc->wavenum;

    if (center == srcRWG->center) return 0.0; // TODO: Handle self-interactions
 
    cmplx intRad = 0.0;

    int obsTriIdx = 0;
    for (const auto& obsTri : tris) {

        const auto& obsQuads = obsTri->quads;
        const auto& obsXnc = Xnc[obsTriIdx];

        int srcTriIdx = 0;
        for (const auto& srcTri : srcRWG->tris) {

            const auto& srcQuads = srcTri->quads;
            const auto& srcXnc = srcRWG->Xnc[srcTriIdx];
            
            if (obsTri == srcTri) {
                ++srcTriIdx;
                continue; // TODO: Handle coincident tris
            }

            for (const auto& [obs, obsWeight] : obsQuads) {
                for (const auto& [src, srcWeight] : srcQuads) {

                    const vec3cd& rad = 
                        srcWeight 
                        * Math::dyadicG(obs-src, k) * (src-srcXnc) 
                        * Math::sign(srcTriIdx);

                    intRad += 
                        obsWeight 
                        * conj(rad.dot(obs-obsXnc)) // Hermitian dot!
                        * Math::sign(obsTriIdx);
                }
            }

            ++srcTriIdx;
        }

        ++obsTriIdx;
    }

    return leng * srcRWG->leng * intRad;
}