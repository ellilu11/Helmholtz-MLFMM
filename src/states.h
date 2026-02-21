#pragma once

#include "types.h"

struct States {
    States() = default;

    States(int nsrcs) 
    {
        lvec = vecXcd::Zero(nsrcs);
        rvec = vecXcd::Zero(nsrcs);
        currents = vecXcd::Zero(nsrcs);
    }

    vecXcd lvec;
    vecXcd rvec;
    vecXcd currents;
};