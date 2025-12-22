#pragma once

#include "math.h"

constexpr double c0 = 299792458.0;
constexpr double mu0 = 1.256637E-6;

const cmplx C = -iu * c0 * mu0 / (4.0 * PI);