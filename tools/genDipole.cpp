#include <filesystem>
#include "../src/dipole.h"

std::filesystem::path makePath(const Config& config) {
    std::string distStr =
        [&]() -> std::string {
        switch (config.pdist) {
            case Dist::UNIFORM:    return "uniform";
            case Dist::GAUSSIAN:   return "gauss";
            case Dist::SPHERE:     return "sphere";
            case Dist::CYLINDER:   return "cyl";
        }
        }();

    return
        std::filesystem::path("config") / "dipole" /
        (distStr + "_n" + std::to_string(config.nsrcs) + ".txt");
}

template <class dist0, class dist1 = dist0, class dist2 = dist0>
void genDipoles(const Config& config, const shared_ptr<Exct::PlaneWave> Einc)
{
    using namespace std;

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

    auto fpath = makePath(config);

    ofstream srcFile(fpath);
    for (const auto& src : srcs) srcFile << dipoles);
}