#pragma once

#include "rwg.h"
#include "subrwg.h"

class Mesh::SrcRWG final : public RWG {
    friend class BC;

public:
    SrcRWG(std::shared_ptr<Excitation::PlaneWave>,
           size_t,
           const Eigen::Vector4i&);

    void findSubRWGs();

    void propagateRvals();

    void buildBC() { bc = std::make_unique<BC>(this); }

    // vec3d getCenter() const override { return glVerts[iCenter]; } // needs mesh to be refined

private:
    std::array<int,10> iSubs;
    std::unique_ptr<BC> bc;
};