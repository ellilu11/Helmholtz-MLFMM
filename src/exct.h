#pragma once

extern const Config config;

namespace Exct {

    struct PlaneWave {
        PlaneWave()
            : pol(vec3d{ 1,0,0 }), 
            wavehat(vec3d{ 0,0,1 }), 
            wavevec(vec3d{ 0,0,1 }), 
            amplitude(1.0)
        {};

        PlaneWave(const vec3d& pol, const vec3d& wavehat, double amplitude)
            : pol(pol.normalized()), 
              wavehat(wavehat.normalized()),
              wavevec(config.k*this->wavehat), 
              amplitude(amplitude)
        {
            std::cout << " Importing excitation...\n";

            if (std::abs(this->pol.dot(this->wavehat)) > 1e-6)
                throw std::runtime_error("Polarization and wave vector must be orthogonal");

            std::cout << "   Polarization:    " << this->pol << '\n';
            std::cout << "   Wave vector:     " << this->wavevec << "\n\n";
        };

        vec3d pol;          // unit polarization
        vec3d wavehat;      // unit wavevector
        vec3d wavevec;      // wavevector
        double amplitude;   // amplitude
    };

    using PlaneWaves = std::vector<std::unique_ptr<PlaneWave>>;
    PlaneWaves Eincs; // incident plane waves
}

