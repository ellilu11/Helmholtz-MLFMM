#pragma once

#include "rwg.h"
#include "subrwg.h"

class SrcRWG final : public RWG {
    friend class BC;

public:
    static SrcVec importRWG(
        const std::filesystem::path&,
        const std::filesystem::path&,
        const std::filesystem::path&,
        const Precision,
        const std::shared_ptr<Excitation::PlaneWave>);

    SrcRWG(std::shared_ptr<Excitation::PlaneWave>,
           size_t,
           const Eigen::Vector4i&,
           const TriVec&);

    void buildSubRWGs();

    void buildBC() { bc = std::make_unique<BC>(this); }

private:
    std::array<std::shared_ptr<SubRWG>,14> subrwgs;
    std::unique_ptr<BC> bc;

    vec2i glIdxc;  // global indices of common vertices
    vec2i glIdxnc; // global indices of non-common vertices
};