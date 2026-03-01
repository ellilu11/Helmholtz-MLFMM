#pragma once

#include <random>
#include "source.h"

class Dipole final : public Source {
public:
    Dipole() = default;

    Dipole(
        std::shared_ptr<Exc::PlaneWave> Einc, size_t iSrc, const vec3d& X)
        : Source(std::move(Einc), iSrc), pos(X), 
        pmag(Phys::p0), pol(vec3d(pmag, 0, 0)), phat(pol/pmag)
    {};

    Dipole(
        std::shared_ptr<Exc::PlaneWave> Einc, size_t iSrc, const vec3d& X, const vec3d& P)
        : Dipole(std::move(Einc), iSrc, X)
    {
        pol = P;
        pmag = P.norm();
        phat = P / pmag;

        buildVoltage();
    };

    void buildVoltage() override {
        voltage = -Einc->amplitude * exp(iu*Einc->wavevec.dot(pos))
            * pol.dot(Einc->pol);
    }

    vec3d getCenter() const override { return pos; }

    /* getRadAlongDir(X,kvec)
     * Return the outgoing radiated amplitude at X along direction kvec
     * due to this dipole
     * X    : observation point (Cartesian)
     * kvec : wavevector
     */
    vec3cd getRadAlongDir(
        const vec3d& X, const vec3d& kvec) const override 
    {
        return exp(iu*kvec.dot(X-pos)) * phat;
    }

    vec3cd getFarAlongDir(
        const vec3d& krhat) const override
    {
        return exp(-iu*pos.dot(krhat)) * phat;
    }

    /* getIntegratedRad(src)
     * Return the radiated field due to src tested with this dipole
     */
    cmplx getIntegratedRad(const std::shared_ptr<Source> src) const override {
        auto srcDip = dynamic_pointer_cast<Dipole>(src);

        if (pos == srcDip->pos) return 0.0; // TODO: Radiation reaction field

        vec3cd rad = Math::dyadicG(pos - srcDip->pos, k) * srcDip->phat;

        return conj(rad.dot(phat));
    }

    friend std::ostream& operator<<(std::ostream& os, Dipole& src) {
        os << src.pos << ' ' << src.pol << '\n';
        return os;
    }

private:
    vec3d pos;   // position
    vec3d pol;   // pol. density vector
    vec3d phat;  // unit pol. density vector
    double pmag; // pol. density magnitude 
};