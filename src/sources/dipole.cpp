#include "dipole.h"

SrcVec importDipoles(
    const std::filesystem::path& fpath, 
    const std::shared_ptr<PlaneWave>& Einc) 
{
    std::ifstream inFile(fpath);
    if (!inFile) throw std::runtime_error("Unable to find file");
    std::string line;
    SrcVec dipoles;

    while (getline(inFile, line)) {
        std::istringstream iss(line);

        vec3d pos;
        if (iss >> pos)
            dipoles.emplace_back(make_shared<Dipole>(Einc, pos));
        else
            throw std::runtime_error("Unable to parse line");
    }
    return dipoles;
}

template <class dist0, class dist1 = dist0, class dist2 = dist0>
SrcVec makeDipoles(const Config& config, const std::shared_ptr<PlaneWave> Einc)
{
    using namespace std;

    SrcVec dipoles;

    random_device rd;
    mt19937 gen(rd());

    dist0 rand0(0, 1);
    dist1 rand1(0, 1);
    dist2 rand2(0, 1);

    if (config.pdist == Dist::UNIFORM) {
        double lim = config.rootLeng/2.0;
        rand0 = dist0(-lim, lim);
        rand1 = dist1(-lim, lim);
        rand2 = dist2(-lim, lim);
    }

    for (int n = 0; n < config.nsrcs; ++n) {
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

        dipoles.push_back(make_shared<Dipole>(Einc, X));
    }

    return dipoles;
}

void Dipole::buildRHS() {

    const auto& kvec = Einc->wavenum * Einc->wavevec;

    rhs = -Einc->amplitude * exp(iu*kvec.dot(pos)) 
            * pol * phat.dot(Einc->pol);

}

void Dipole::buildCurrent() {

    current = iu * c0 * Einc->wavenum * pol; // |J| = i \omega |P|

    // TODO: Predict current from rhs vector
}

/* getRadAlongDir(X,kvec)
 * Return the outgoing radiated amplitude due to this
 * dipole at X along direction kvec
 * X    : observation point (Cartesian)
 * kvec : wavevector
 */
vec3cd Dipole::getRadAlongDir(
    const vec3d& X, const vec3d& kvec) const {

    return current * exp(iu*kvec.dot(X-pos)) * phat;
}

/* getIncAlongDir(X,kvec)
 * Return the incoming radiated amplitude at this
 * dipole due to field at X along direction kvec
 * X    : source point (Cartesian)
 * kvec : wavevector
 */
vec3cd Dipole::getIncAlongDir(
    const vec3d& X, const vec3d& kvec) const {

    return current * exp(iu*kvec.dot(pos-X)) * phat;
}

/* getRad(X,k)
 * Return the full radiated field with given wavenum
 * due to this dipole at X
 */
vec3cd Dipole::getRad(const vec3d& X, double wavenum) const {

    const auto& dyadic = Math::dyadicG(X - pos, wavenum);

    return current * dyadic * phat;
}

/* getIntegratedRad(X,k)
 * Return the full radiated field with given wavenum
 * due to src tested at this dipole
 */
cmplx Dipole::getIntegratedRad(const std::shared_ptr<Source> src, double wavenum) const {

    auto srcDip = dynamic_pointer_cast<Dipole>(src);

    return srcDip->getRad(pos, wavenum).dot(phat);
}



