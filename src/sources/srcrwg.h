#pragma once

#include "rwg.h"
#include "subrwg.h"

class SrcRWG final : public RWG {
    friend class BC;

public:
    SrcRWG(std::shared_ptr<Excitation::PlaneWave>,
           size_t,
           const Eigen::Vector4i&);

    void buildSubRWGs();

    void buildBC() { bc = std::make_unique<BC>(this); }

private:
    std::array<std::shared_ptr<SubRWG>,14> subrwgs;
    std::unique_ptr<BC> bc;
};