#pragma once

#include "fmm.h"

struct FMM::Coeffs {
    std::vector<cmplx> vals;
    size_t nDir;

    Coeffs() = default;

    Coeffs(size_t nDir, cmplx val = 0.0)
        : vals(2*nDir, val), nDir(nDir)
    {}

    size_t size() const { return nDir; }

    void resize(size_t nDir) {
        vals.resize(2.0*nDir);
        this->nDir = nDir;
    }

    void fill(cmplx val) {
        std::fill(vals.begin(), vals.end(), val);
    }

    void setCoeffAlongDir(const vec2cd& vec, size_t iDir) {
        vals[iDir] = vec[0];
        vals[nDir+iDir] = vec[1];
    }

    vec2cd getVecAlongDir(size_t iDir) const {
        return vec2cd(vals[iDir], vals[nDir+iDir]);
    }

    Coeffs& operator+=(const Coeffs& coeffs) noexcept {
        vals += coeffs.vals;
        return *this;
    }

    Coeffs& operator*=(cmplx val) noexcept {
        vals *= val;
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

    friend std::ostream& operator<<(std::ostream& os, const Coeffs& coeffs) noexcept {
        os << coeffs.vals;
        return os;  
    }
};