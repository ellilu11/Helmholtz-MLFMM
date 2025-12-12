#pragma once

#include "math.h"

struct Source {
    Source()
        : 
        amplitude(1.0), 
        wavenum(1.0),
        pol(vec3d{ 1,0,0 }),
        wavevec(vec3d{ 0,0,1 })
    {
    };

    /* TODO: Read src params from file
    Src(const std::string& fileName) {
        std::ifstream is(fileName);
        is >> amplitude >> wavenum >> pol >> wavevec;
    }*/

    double amplitude;   // amplitude
    double wavenum;     // wavenumber
    vec3d pol;          // unit polarization
    vec3d wavevec;      // unit wavevector
};