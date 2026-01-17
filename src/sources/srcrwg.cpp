#include "srcrwg.h"

SrcVec SrcRWG::importRWG(
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

    while (std::getline(file,line)) {
        std::istringstream iss(line);
        Eigen::Vector4i idx4;

        if (iss >> idx4)
            rwgs.push_back(std::make_shared<SrcRWG>(Einc,rwgIdx++,idx4,triangles));
        else
            throw std::runtime_error("Unable to parse line");
    }

    SubRWG::buildVertsToSubrwgs(Triangle::glVerts.size());

    for (const auto& rwg : rwgs)
        dynamic_pointer_cast<SrcRWG>(rwg)->buildBC();

    return rwgs;
}

SrcRWG::SrcRWG(
    std::shared_ptr<Excitation::PlaneWave> Einc,
    size_t rwgIdx,
    const Eigen::Vector4i& idx4,
    const TriVec& triangles)
    : RWG(std::move(Einc), rwgIdx, std::move(triangles[idx4[2]]), std::move(triangles[idx4[3]])),
    // tris({ triangles[idx4[2]], triangles[idx4[3]] }),
    glIdxc({ idx4[0],idx4[1] })
{
    Xc[0] = Triangle::glVerts[glIdxc[0]], Xc[1] = Triangle::glVerts[glIdxc[1]];
    center = (Xc[0]+Xc[1])/2.0;
    leng = (Xc[0]-Xc[1]).norm();

    // Find non-common vertices
    for (int i = 0; i < 2; ++i)
        for (const auto& glIdxs : tris[i]->glIdxs)
            if (glIdxs != glIdxc[0] && glIdxs != glIdxc[1]) {
                glIdxnc[i] = glIdxs;
                Xnc[i] = Triangle::glVerts[glIdxs];
            }

    buildSubRWGs();

    buildVoltage(); // needs Xnc initialized!

    //std::cout << glIdxc[0] << ' '<< glIdxc[1] << ' '
    //    << glIdxc_pm[0] << ' ' << glIdxc_pm[1] << "\n";

    //std::cout << '(' << X[0] << ") (" << X[1] << ") ("
    //    << Xnc[0] << ") (" << Xnc[1] << ") " << leng << '\n';
};

void SrcRWG::buildSubRWGs() {
    const Eigen::MatrixXi subIdxs{
        {0, 0, 4, 1, 2, 3},
        {5, 1, 5, 2, 3, 4}
    };

    int iTri = 0;
    TriVec midSubtris;

    for (const auto& tri : tris) {
        const auto& subtris = tri->getSubtris(vec3i(glIdxnc[iTri],glIdxc[0],glIdxc[1]));
        const int iTri8 = 8*iTri;

        for (int iSub = 0; iSub < 6; ++iSub) {
            const auto& subtris0 = subtris[subIdxs(0,iSub)];
            const auto& subtris1 = subtris[subIdxs(1,iSub)];

            auto subrwg = std::make_shared<SubRWG>(subtris0,subtris1);

            subrwgs[iTri8+iSub] = subrwg; // TODO: Move assign
            glSubrwgs.push_back(std::move(subrwg));

            if (iSub == 2 || iSub == 3) midSubtris.push_back(subtris[iSub]);
        }

        ++iTri;
    }

    // Construct RWGs along common edge of parent RWG
    for (int iSub = 0; iSub < 2; ++iSub) {
        const auto& subtris0 = midSubtris[iSub];
        const auto& subtris1 = midSubtris[iSub+2];

        auto subrwg = std::make_shared<SubRWG>(subtris0,subtris1);
        subrwgs[iSub+6] = subrwg; // TODO: Move assign
        glSubrwgs.push_back(std::move(subrwg));
    }

    //for (const auto& rwg : subrwgs) {
    //    for (auto& idx : rwg->glIdxs)
    //        std::cout << idx << ' ';
    //    std::cout << '\n';
    //}
        // std::cout << rwg->Xc[0] << ' ' << rwg->Xc[1] << '\n';
}