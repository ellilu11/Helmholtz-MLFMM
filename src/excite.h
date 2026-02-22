#pragma once

#include "types.h"

extern double k;

namespace Exc {
    struct PlaneWave {
        PlaneWave()
            : pol(vec3d{ 1,0,0 }),
            wavevec(vec3d{ 0,0,1 }),
            wavenum(wavevec.norm()),
            amplitude(1.0)
        {};

        PlaneWave(
            const vec3d& pol,
            const vec3d& wavehat,
            double wavenum, double amplitude)
            : pol(pol),
            wavevec(wavenum*wavehat),
            wavenum(wavenum),
            amplitude(amplitude)
        {};

        vec3d pol;          // unit polarization
        vec3d wavevec;      // wavevector
        double wavenum;     // wavenumber
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
    if (iss >> amplitude >> pol >> wavehat >> wavenum)
        Einc = std::make_shared<Exc::PlaneWave>(pol, wavehat, wavenum, amplitude);
    else
        throw std::runtime_error("Unable to parse line");

    ::k = Einc->wavenum; // set global wavenumber
    std::cout << "   Polarization:    " << pol << '\n';
    std::cout << "   Unit wavenum:    " << wavehat << '\n';
    std::cout << "   Wave number:     " << k << "/m\n\n";

    return Einc;
}

