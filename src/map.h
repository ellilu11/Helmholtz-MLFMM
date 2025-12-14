#pragma once

#include <map>
#include <unordered_map>

constexpr double DistEPS = 1.0E-3;
constexpr double PsiEPS = 1.0E-9;

struct Comp {
    bool operator()(double x, double y) const noexcept {
        return x + DistEPS < y;
    }
};

struct HashFunc {
    std::size_t operator()(double x) const noexcept {
        long long q = static_cast<long long>(std::llround(x / PsiEPS));
        return std::hash<long long>{}(q);
    }
};

struct HashEq {
    bool operator()(double x, double y) const noexcept {
        return std::fabs(x - y) <= PsiEPS;
    }
};

template <typename T>
using Map = std::map<double, T, Comp>;

template <typename T>
using HashMap = std::unordered_map<double, T, HashFunc, HashEq>;