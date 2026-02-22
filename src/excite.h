#pragma once

#include "types.h"

extern double k;

namespace Exc {
    struct PlaneWave {
        PlaneWave()
            : pol(vec3d{1,0,0}), wavevec(vec3d{0,0,1}), amplitude(1.0)
        {};

        PlaneWave(const vec3d& pol, const vec3d& wavehat, double amplitude)
            : pol(pol.normalized()), wavevec(k*wavehat.normalized()), amplitude(amplitude)
        {
            std::cout << "   Polarization:    " << this->pol << '\n';
            std::cout << "   Wave vector:     " << this->wavevec << '\n';
            std::cout << "   Amplitude:       " << this->amplitude << "\n\n";
        };

        vec3d pol;          // unit polarization
        vec3d wavevec;      // wavevector
        double amplitude;   // amplitude
    };

    std::shared_ptr<PlaneWave> importPlaneWaves(const std::filesystem::path&);
}

std::shared_ptr<Exc::PlaneWave> 
    Exc::importPlaneWaves(const std::filesystem::path& fpath)
{
    std::ifstream inFile(fpath);
    if (!inFile) throw std::runtime_error("Unable to find file");

    std::string line;
    std::getline(inFile, line);
    std::istringstream iss(line);

    std::shared_ptr<Exc::PlaneWave> Einc;

    vec3d pol, wavehat;
    double wavenum, amplitude;
    if (iss >> pol >> wavehat >> amplitude)
        Einc = std::make_shared<Exc::PlaneWave>(pol, wavehat, amplitude);
    else
        throw std::runtime_error("Unable to parse line");

    return Einc;
}

