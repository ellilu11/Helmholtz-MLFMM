#include "srcrwg.h"

SrcRWG::SrcRWG(
    std::shared_ptr<Excitation::PlaneWave> Einc,
    size_t iSrc,
    const Eigen::Vector4i& idx4)
    : RWG(std::move(Einc), iSrc, idx4)
{
    // buildSubRWGs();

    buildVoltage(); 

    //std::cout << "Built srcRWG #" << iSrc << " w/ common vertices # " 
    //    << iVertsC[0] << ' '<< iVertsC[1] << " and non-common vertices # "
    //    << iVertsNC[0] << ' ' << iVertsNC[1] << "\n";

    //std::cout << '(' << X[0] << ") (" << X[1] << ") ("
    //    << Xnc[0] << ") (" << Xnc[1] << ") " << leng << '\n';
};

void SrcRWG::buildSubRWGs() {
    const Eigen::MatrixXi subIdxs{
        {0, 0, 4, 1, 2, 3},
        {5, 1, 5, 2, 3, 4}
    };

    int iTri = 0;
    std::vector<Triangle> midSubtris;
    /*
    for (const auto& tri : tris) {
        const auto& subtris = tri->getSubtris(vec3i(iVertsNC[iTri],iVertsC[0],iVertsC[1]));
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
    */

    //for (const auto& rwg : subrwgs) {
    //    for (auto& idx : rwg->glIdxs)
    //        std::cout << idx << ' ';
    //    std::cout << '\n';
    //}
        // std::cout << rwg->Xc[0] << ' ' << rwg->Xc[1] << '\n';
}