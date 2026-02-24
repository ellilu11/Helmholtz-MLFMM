#pragma once

#include <cctype>
#include <fstream>
#include <iostream>
#include <type_traits>
#include "phys.h"

extern double k;

enum class Mode { READ, WRITE };

enum class Precision { VERYLOW, LOW, MEDIUM, HIGH, VERYHIGH };

enum class Dist { UNIFORM, GAUSSIAN, SPHERE, CYLINDER };

enum class QDist { UNIFORM, RANDSIGN, RANDOM };

void getDigit(std::istringstream& iss, char ch) {
    while (iss.get(ch)) {
        if (std::isdigit(static_cast<unsigned char>(ch))) {
            iss.unget();
            break;
        }
    }
}

template <typename T>
std::ifstream& operator>>(std::ifstream& is, T& val) {
    std::string line;
    if (std::getline(is, line)) {
        std::istringstream iss(line);

        char ch = '\0';
        getDigit(iss, ch);

        if constexpr (std::is_enum_v<T>) {
            typename std::underlying_type<T>::type eval;

            while (iss >> eval) {
                val = static_cast<T>(eval);
                getDigit(iss, ch);
            }

        } else
            while (iss >> val)
                getDigit(iss, ch);
    }

    return is;
}

int getNumQuads(Precision prec) {
    return [&]() {
        switch (prec) {
            case Precision::VERYLOW:  return 1;
            case Precision::LOW:      return 3;
            case Precision::MEDIUM:   return 7;
            case Precision::HIGH:     return 13;
            case Precision::VERYHIGH: return 19;
        };
        } ();
}

struct Config {
    Config() = default;
    Config(const std::string& fileName) {
        std::ifstream is(fileName);
        is >> mode >> pdist >> qdist >> quadPrec // enums first
            >> digits >> interpOrder >> overInterp
            >> rootLeng >> maxNodeSrcs  >> evalDirect
            >> nsrcs >> wavenum;

        ::k = wavenum; // set global wavenumber

        std::cout << " *********************** \n";
        std::cout << " *** Helmholtz-MLFMM *** \n";
        std::cout << " *********************** \n";

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "   Mode:            " << (mode == Mode::READ ? "READ" : "WRITE") << '\n';
        std::cout << "   # Sources:       " << nsrcs << '\n';
        std::cout << "   Digit precision: " << digits << '\n';
        std::cout << "   Interp order:    " << interpOrder << '\n';
        std::cout << "   Overinterp:      " << overInterp << '\n';
        std::cout << "   Max node srcs:   " << maxNodeSrcs << '\n';
        std::cout << "   Tri quad rule:   " << getNumQuads(quadPrec) << "-point\n";
        std::cout << "   Root length:     " << rootLeng << " m\n";
        std::cout << "   Wavenumber:      " << k << " /m\n\n";
    }

    Mode mode;
    Precision quadPrec;
    int digits;
    int interpOrder;
    double overInterp;
    double rootLeng;
    int maxNodeSrcs;
    bool evalDirect;
    int nsrcs;
    double wavenum;

    // For point dipoles only
    Dist pdist;
    QDist qdist;
};