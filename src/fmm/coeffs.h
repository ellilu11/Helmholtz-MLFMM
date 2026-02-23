#pragma once

#include "fmm.h"

struct FMM::Coeffs {
    std::vector<cmplx> theta;
    std::vector<cmplx> phi;

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

    void setCoeffAlongDir(const vec2cd& vec, size_t iDir) {
        theta[iDir] = vec[0];
        phi[iDir] = vec[1];
    }

    vec2cd getVecAlongDir(size_t iDir) const {
        return vec2cd(theta[iDir], phi[iDir]);
    }

    /*
    std::vector<cmplx> flattened() const {
        std::vector<cmplx> flatVec = theta;
        flatVec.insert(flatVec.end(), phi.begin(), phi.end());

        return flatVec;
    }*/

    size_t size() const { return theta.size(); }

    Coeffs& operator+=(const Coeffs& coeffs) noexcept {
        theta += coeffs.theta;
        phi += coeffs.phi;
        return *this;
    }

    Coeffs& operator*=(cmplx val) noexcept {
        theta *= val;
        phi *= val;
        return *this;
    }

    friend Coeffs operator+(Coeffs coeffs0, const Coeffs& coeffs1) noexcept {
        coeffs0 += coeffs1;
        return coeffs0;
    }

    friend Coeffs operator*(cmplx val, Coeffs coeffs) noexcept {
        coeffs *= val;
        return coeffs;
    }

    friend std::ostream& operator<<(std::ostream& os, const Coeffs& coeffs) {
        os << coeffs.theta;
        os << coeffs.phi;
        return os;  
    }
};