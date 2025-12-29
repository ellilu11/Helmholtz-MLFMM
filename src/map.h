#pragma once

#include <map>
#include <unordered_map>

namespace Maps {

    // Pick largest values that avoid collisions
    constexpr double DistEPS = 1.0E-2;
    constexpr double EPS = 1.0E-3; 
    constexpr double VEPS = 3.0*EPS;

    struct DoubleComp {
        bool operator()(double x, double y) const noexcept {
            return x + DistEPS < y;
        }
    };

    struct DoubleFunc {
        std::size_t operator()(double x) const noexcept {
            long long q = static_cast<long long>(std::llround(x / EPS));
            return std::hash<long long>{}(q);
        }
    };

    struct DoubleEq {
        bool operator()(double x, double y) const noexcept {
            return std::fabs(x - y) <= EPS;
        }
    };

    struct VecFunc {
        std::size_t operator()(const vec3d& X) const noexcept {
            long long q = static_cast<long long>(std::llround(X.norm() / VEPS));
            return std::hash<long long>{}(q);
        }
    };

    struct VecEq {
        bool operator()(const vec3d& X, const vec3d& Y) const noexcept {
            return ((X-Y).norm()) <= VEPS;
        }
    };
}

// TODO: Concepts
template <typename T>
using Map = std::map<double, T, Maps::DoubleComp>;

template <typename T>
using HashMap = std::unordered_map<double, T, Maps::DoubleFunc, Maps::DoubleEq>;

template <typename T>
using VecHashMap = std::unordered_map<vec3d, T, Maps::VecFunc, Maps::VecEq>;