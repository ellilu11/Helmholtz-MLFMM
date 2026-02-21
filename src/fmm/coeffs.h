#pragma once

#include "types.h"

struct FMM::Coeffs {
    cmplxVec theta;
    cmplxVec phi;

    Coeffs() = default;

    Coeffs(size_t N, cmplx val = 0.0) 
        : theta(N, val), phi(N, val) 
    {}

    void resize(size_t N) {
        theta.resize(N, 0.0);
        phi.resize(N, 0.0);
    }

    void fillZero() {
        std::fill(theta.begin(), theta.end(), 0.0);
        std::fill(phi.begin(), phi.end(), 0.0);
    }

    vec2cd getVecAlongDir(size_t idx) const {
        return vec2cd(theta[idx], phi[idx]);
    }

    size_t size() const { return theta.size(); }

    Coeffs& operator+=(const Coeffs& coeffs) noexcept {
        theta += coeffs.theta;
        phi += coeffs.phi;
        return *this;
    }

    friend Coeffs operator+(Coeffs coeffs0, const Coeffs& coeffs1) noexcept {
        coeffs0 += coeffs1;
        return coeffs0;
    }
};