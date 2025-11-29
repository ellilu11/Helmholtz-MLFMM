#pragma once

#include "math.h"

struct Src {
    Src()
        : amplitude(1.0), k(1.0),
        pol(vec3d{ 1,0,0 }),
        wavevec(vec3d{ 0,0,1 })
    {
    };

    double amplitude;   // amplitude
    double k;           // wavenumber
    vec3d pol;          // unit polarization
    vec3d wavevec;      // unit wavevector
};