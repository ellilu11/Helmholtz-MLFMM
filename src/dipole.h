#pragma once

#include <random>
#include "source.h"

class Dipole final : public Source {
public:
    Dipole() = default;

    Dipole(const vec3d& X, size_t iSrc)
        : Source(iSrc), pos(X), 
        pmag(Phys::p0), pol(vec3d(pmag, 0, 0)), phat(pol/pmag)
    {};

    Dipole(const vec3d& X, const vec3d& P, size_t iSrc)
        : Dipole(iSrc, X)
    {
        pol = P;
        pmag = P.norm();
        phat = P / pmag;
    };

    cmplx getVoltage() override {
        using namespace Exct;

        cmplx voltage = 0.0;
        for (const auto& Einc : Eincs)
            voltage -= Einc->amplitude * exp(iu*Einc->wavevec.dot(pos))
                * pol.dot(Einc->pol);

        return voltage;
    }

    vec3d getCenter() const override { return pos; }

    std::pair<vec3cd, vec3cd> getRadsAlongDir(const vec3d& X, const vec3d& kvec) const override
    {
        vec3cd rad = exp(iu*kvec.dot(X-pos)) * phat / (4.0*PI); // apply factor of 1/(4pi)
        return { rad, kvec.cross(rad).conjugate() }; // Hermitian cross!
    }

    vec3cd getFarAlongDir(
        const vec3d& krhat) const override
    {
        return exp(-iu*pos.dot(krhat)) * phat / (4.0*PI); // apply factor of 1/(4pi)
    }

    /* getIntegratedEFIE(src)
     * Return the electric field due to src tested with this dipole
     */
    cmplx getIntegratedEFIE(const std::shared_ptr<Source> src, double k) const override {
        auto srcDip = dynamic_pointer_cast<Dipole>(src);
        if (pos == srcDip->pos || config.alpha == 0.0) return 0.0;

        vec3cd rad = Math::dyadicG(pos - srcDip->pos, k) * srcDip->phat;

        return phat.dot(rad);
    }

    /* getIntegratedMFIE(src)
     * Return the magnetic field due to src tested with this dipole
     */
    cmplx getIntegratedMFIE(const std::shared_ptr<Source> src, double k) const override {
        auto srcDip = dynamic_pointer_cast<Dipole>(src);
        if (pos == srcDip->pos || config.alpha == 1.0) return 0.0;
        vec3d R = pos - srcDip->pos;
        double r = R.norm(), r3 = r*r*r;

        vec3cd gradG = (-1.0+iu*k*r)*exp(iu*k*r)/r3 * R;
        vec3cd rad = -(gradG.cross(srcDip->phat)).conjugate(); // Hermitian cross!

        return phat.dot(rad);
    }

    double getIntegratedMass(const std::shared_ptr<Source> src) const override {
        return 0.0; // Dipole magnetic self interaction?
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