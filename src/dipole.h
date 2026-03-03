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
        vec3d polInc =
            config.alpha * Einc->pol +
            config.beta * Einc->wavehat.cross(Einc->pol);

        voltage = -Einc->amplitude * exp(iu*Einc->wavevec.dot(pos))
            * pol.dot(polInc);
    }

    vec3cd getRadAlongDir(const vec3d& X, const vec3d& kvec) const override {
        //vec3cd intPlaneWave = exp(iu*kvec.dot(X-pos)) * phat;
        //return config.alpha * intPlaneWave -
        //    config.beta * iu * kvec.cross(intPlaneWave);
        return exp(iu*kvec.dot(X-pos)) * phat;
    }

    vec3cd getFarAlongDir(const vec3d& krhat) const override {
        return exp(-iu*pos.dot(krhat)) * phat;
    }

    /* getIntegratedEFIE(src)
     * Return the electric field due to src tested with this dipole
     */
    cmplx getIntegratedEFIE(const std::shared_ptr<Source> src) const override {
        auto srcDip = dynamic_pointer_cast<Dipole>(src);
        if (pos == srcDip->pos || config.alpha == 0.0) return 0.0;

        vec3cd rad = Math::dyadicG(pos - srcDip->pos, k) * srcDip->phat;
        
        return phat.dot(rad);
    }

    /* getIntegratedMFIE(src)
     * Return the magnetic field due to src tested with this dipole
     */
    cmplx getIntegratedMFIE(const std::shared_ptr<Source> src) const override {
        auto srcDip = dynamic_pointer_cast<Dipole>(src);
        if (pos == srcDip->pos || config.beta == 0.0) return 0.0;
        vec3d R = pos - srcDip->pos;
        double r = R.norm(), r3 = r*r*r;

        vec3cd gradG = (-1.0+iu*k*r)*exp(iu*k*r)/r3 * R;
        vec3cd rad = -gradG.cross(srcDip->phat);

        return phat.dot(rad);
    }

    vec3d getCenter() const override { return pos; }

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