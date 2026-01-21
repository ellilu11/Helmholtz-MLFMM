#pragma once

#include "rwg.h"
#include "subrwg.h"

class SrcRWG final : public RWG {
    friend class BC;

public:
    SrcRWG(std::shared_ptr<Excitation::PlaneWave>,
           size_t,
           const Eigen::Vector4i&);

    void buildSubIdx();

    void buildBC() { bc = std::make_unique<BC>(this); }

    vec3d getCenter() const override { return Triangle::glVerts[iCenter]; }

private:
    std::array<int,10> iSubs;
    std::unique_ptr<BC> bc;

    int iCenter; // index of midpoint of common edge
};