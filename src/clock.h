#pragma once

#include <chrono>

using Time = std::chrono::duration<double, std::milli>;
using Clock = std::chrono::high_resolution_clock;

struct ClockTimes {
    ClockTimes() = default;

    Time S2M{};
    Time M2M{};
    Time M2L{};
    Time L2L{};
    Time L2T{};
    Time S2T{};
};