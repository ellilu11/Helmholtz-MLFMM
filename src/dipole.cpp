#include "dipole.h"

using namespace std;

template <class dist0, class dist1 = dist0, class dist2 = dist0>
SrcVec makeDipoles(const Config& config, const shared_ptr<Exc::PlaneWave> Einc)
{
    SrcVec dipoles;

    random_device rd;
    mt19937 gen(rd());

    dist0 rand0(0, 1);
    dist1 rand1(0, 1);
    dist2 rand2(0, 1);

    dist0 prand0(0, 1);
    dist1 prand1(0, 1);

    if (config.pdist == Dist::UNIFORM) {
        double lim = config.rootLeng/2.0;
        rand0 = dist0(-lim, lim);
        rand1 = dist1(-lim, lim);
        rand2 = dist2(-lim, lim);
    }

    for (size_t n = 0; n < config.nsrcs; ++n) {
        vec3d X = [&] {
            double r, th, ph, z;

            switch (config.pdist) {
                case Dist::UNIFORM:
                    return vec3d(rand0(gen), rand1(gen), rand2(gen));

                    // 2D Gaussian at z = 0; TODO: 3D Gaussian
                case Dist::GAUSSIAN: {
                    r = sqrt(-2.0 * log(rand0(gen)));
                    th = PI / 2.0;
                    ph = 2.0 * PI * rand1(gen);
                    return Math::fromSph(vec3d(r, th, ph));
                }

                case Dist::SPHERE: {
                    r = 0.90 * (config.rootLeng / 2.0);
                    th = acos(2.0*rand0(gen) - 1.0);
                    ph = 2.0 * PI * rand1(gen);
                    return Math::fromSph(vec3d(r, th, ph));
                }

                case Dist::CYLINDER: {
                    r = 0.45 * (config.rootLeng / 2.0);
                    ph = 2.0 * PI * rand1(gen);
                    z = 0.90 * config.rootLeng * (rand1(gen) - 0.5);
                    return Math::fromCyl(vec3d(r, ph, z));
                }
            }
            }();

        vec3d P = [&] {
            double th, ph;

            switch (config.qdist) {
                case QDist::UNIFORM:
                    return vec3d(Phys::p0, 0.0, 0.0);

                case QDist::RANDSIGN: {
                    uniform_int_distribution randi(0, 1);
                    return vec3d(Math::sign(randi(gen))*Phys::p0, 0.0, 0.0);
                }

                case QDist::RANDOM: {
                    th = acos(2.0*prand0(gen) - 1.0);
                    ph = 2.0 * PI * prand1(gen);
                    return Math::fromSph(vec3d(Phys::p0, th, ph));
                }
            }
            }();

        dipoles.emplace_back(make_shared<Dipole>(Einc, n, X, P));
    }

    return dipoles;
}

// TODO: Make Dipole static method
SrcVec importDipoles(
    const filesystem::path& fpath,
    const shared_ptr<Exc::PlaneWave>& Einc)
{
    ifstream inFile(fpath);
    if (!inFile) throw runtime_error("Unable to find file");
    string line;
    SrcVec dipoles;
    size_t idx = 0;

    while (getline(inFile, line)) {
        istringstream iss(line);

        vec3d pos, dip;
        if (iss >> pos >> dip)
            dipoles.emplace_back(make_shared<Dipole>(Einc, idx++, pos, dip));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return dipoles;
}