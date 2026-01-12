#include "bc.h"

BC::BC(std::shared_ptr<Excitation::PlaneWave> Einc,
    size_t idx,
    vec3d X0, vec3d X) :
    Source(std::move(Einc), idx),
    X0(X0), X1(X1)
{
    // Find all refined tris with X0 as vertex and add to member tris

};