#pragma once

#include "fmm.h"

struct FMM::Angles {
    Angles() = default;

    Angles(double, double, int, int);

    void printAngles(std::ofstream&, std::ofstream&);

    pair2i getNumAngles() const {
        return std::make_pair(thetas.size(), phis.size());
    }

    size_t getNumAllAngles() const {
        return thetas.size() * phis.size();
    }

    realVec thetas;       // theta samples
    realVec thetaWeights; // weights of theta samples
    realVec phis;         // phi samples
    int L;                // M2L series truncation number
};

FMM::Angles::Angles(double wavenum, double rootLeng, int digits, int level)
{
    // const double wavenum = Node::wavenum;
    const double nodeLeng = rootLeng / pow(2.0, level);

    // Use excess bandwidth formula
    const int tau = ceil((1.73*wavenum*nodeLeng +
        2.16*pow(digits, 2.0/3.0)*pow(wavenum*nodeLeng, 1.0/3.0)));

    L = floor(0.50*tau); // TODO: Find optimal formula

    // Construct thetas
    const int nth = tau+1;
    std::tie(thetas, thetaWeights) = Interp::gaussLegendre(nth, EPS_NR, 0.0, PI);

    // Absorb sin(theta) into weights
    std::transform(thetaWeights.begin(), thetaWeights.end(), thetas.begin(), thetaWeights.begin(),
        [](double weight, double theta) { return weight * sin(theta); }
    );

    // Construct phis
    const int nph = 2*nth;
    phis.resize(nph);

    for (int iph = 0; iph < nph; ++iph)
        phis[iph] = 2.0*PI*iph/static_cast<double>(nph);

    std::cout << "   (" << level << "," << thetas.size() << "," << phis.size() << ")\n";
}

void FMM::Angles::printAngles(std::ofstream& thfile, std::ofstream& phfile) {
    for (const auto& theta : thetas)
        thfile << theta << '\n';

    for (const auto& phi : phis)
        phfile << phi << '\n';
}