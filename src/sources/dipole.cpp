#include "dipole.h"

Dipole::Dipole(
    std::shared_ptr<PlaneWave> Einc, 
    const vec3d& X) // TODO: Pass in pol. density vector
    : Source(Einc), pos(X), phat(vec3d(1, 0, 0)), pol(1.0)
{
    buildCurrent();
};

void Dipole::buildRHS() {

    const auto& kvec = Einc->wavenum * Einc->wavevec;

    rhs = -Einc->amplitude * exp(iu*kvec.dot(pos)) 
            * pol * phat.dot(Einc->pol);

}

void Dipole::buildCurrent() {

    // current = iu * c0 * Einc->wavenum * pol; // |J| = i \omega |P|
    current = pol;


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

    // std::cout << current << ' ' << exp(iu*kvec.dot(X-pos)) << ' ' << phat << '\n';

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

/* getRadAtPoint(X,k)
 * Return the radiated field due to this dipole 
 * at field point X
 */
vec3cd Dipole::getRadAtPoint(const vec3d& X) const {
    assert(X != pos);

    const auto& dyadic = Math::dyadicG(X - pos, Einc->wavenum);

    // std::cout << dyadic * phat << '\n';

    return current * dyadic * phat;
}

/* getIntegratedRad(X,k)
 * Return the radiated field with given wavenum
 * due to src tested at this dipole
 */
cmplx Dipole::getIntegratedRad(const std::shared_ptr<Source> src) const {

    const auto srcDip = dynamic_pointer_cast<Dipole>(src);

    return srcDip->getRadAtPoint(pos).dot(phat);
}



