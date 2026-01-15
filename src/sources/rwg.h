#pragma once

#include "source.h"
#include "triangle.h"

class RWG;

using RWGVec = std::vector<std::shared_ptr<RWG>>;

class RWG final : public Source {
public:
    friend class BC;

    static SrcVec importRWG(
        const std::filesystem::path&,
        const std::filesystem::path&,
        const std::filesystem::path&,
        const Precision,
        const std::shared_ptr<Excitation::PlaneWave>);

    RWG(std::shared_ptr<Excitation::PlaneWave>,
        size_t,
        const Eigen::Vector4i&, 
        const TriVec&);

    RWG(std::shared_ptr<Triangle>,
        std::shared_ptr<Triangle>);

    void buildSubRWGs();

    static void buildVertsToSubrwgs(int numVerts);

    vec3cd getIntegratedPlaneWave(const vec3d&,bool = 0) const;

    cmplx getIntegratedRad(const std::shared_ptr<Source>) const override;

    double getLeng() const { return leng; }

    vec3d getCenter() const override { return center; } 

    void buildVoltage() override {
        voltage = -Einc->amplitude
            * conj(getIntegratedPlaneWave(Einc->wavevec).dot(Einc->pol)); // Hermitian dot!
    }

    vec3cd getRadAlongDir(const vec3d& X, const vec3d& kvec) const override {
        return exp(iu*kvec.dot(X)) * getIntegratedPlaneWave(kvec).conjugate();
    }

    vec3cd getFarAlongDir(const vec3d& krhat) const override {
        return getIntegratedPlaneWave(krhat).conjugate();
    }

private:
    static RWGVec glSubrwgs; // list of all subrwgs in mesh
    // list of indices of subrwg having vert as common vert
    static std::vector<std::vector<int>> vertsToSubrwgs; 

    std::array<std::shared_ptr<Triangle>,2> tris;

    std::array<vec3d,2> Xc;  // common vertices
    vec2i idx_c; // global indices of common vertices

    std::array<vec3d,2> Xnc;  // non-common vertices
    vec2i idx_nc; // global indices of non-common vertices

    vec3d center;   // midpoint of common edge
    double leng;    // length of common edge

    std::array<std::shared_ptr<RWG>,14> subrwgs;
    std::unique_ptr<BC> bc;
};