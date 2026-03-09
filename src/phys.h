#pragma once

#include "math.h"

namespace Phys {
    constexpr double c0 = 299792458.0;
    constexpr double mu0 = 1.256637E-6;
    constexpr double eta = c0 * mu0; // intrinsic impedance of free space
    constexpr double p0 = 1.0E-8; // Dipole moment
}