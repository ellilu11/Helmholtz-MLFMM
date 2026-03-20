#pragma once

#include "fileio.h"
#include "phys.h"

enum class Mode { FMM, DIR, FMMDIR };

enum class IE { EFIE, MFIE, CFIE };

enum class Precision { VERYLOW, LOW, MEDLOW, MEDIUM, MEDHIGH, HIGH, VERYHIGH};

std::string getModeStr(Mode mode) {
    return [&]() {
        switch (mode) {
            case Mode::FMM:    return "FMM";
            case Mode::DIR:    return "DIRECT";
            case Mode::FMMDIR: return "FMM+DIRECT";
        };
        } ();
}

std::string getIEStr(IE ie) {
    return [&]() {
        switch (ie) {
            case IE::EFIE: return "EFIE";
            case IE::MFIE: return "MFIE";
            case IE::CFIE: return "CFIE";
        };
        } ();
}

size_t getNumQuads(Precision prec) {
    return [&]() {
        switch (prec) {
            case Precision::VERYLOW:  return 1;
            case Precision::LOW:      return 3;
            case Precision::MEDLOW:   return 4;
            case Precision::MEDIUM:   return 6;
            case Precision::MEDHIGH:  return 7;
            case Precision::HIGH:     return 12;
            case Precision::VERYHIGH: return 13;
        };
    } ();
}

struct Config {
    Config() = default;
    
    Config(const std::filesystem::path&);
    
    static std::array<double,14> importConfig(const std::filesystem::path&);

    // General
    Mode mode;          // FMM, direct, or FMM+direct
    Precision quadPrec; // Triangle quadrature precision
    double alpha;       // CFIE parameter
    double k;           // Wavenumber
    // double leps;     // TODO: Length error tolerance for nearfield

    // FMM
    int maxNodeSrcs;    // Max # sources in leaf
    int digits;         // Digit precision
    int interpOrder;    // Interpolation order
    double overInterp;  // Translation operator angular oversampling factor

    // GMRES
    double epsIter;     // Error tolerance for GMRES
    int maxIter;        // Max GMRES iterations (0 for one MVM)
    double iluTol;      // Drop tolerance for ILU preconditioner
    int iluFactor;      // Fill factor for ILU preconditioner

    // CFIE (from alpha)
    IE ie;              // Integral equation type
    cmplx C_efie;       // -eta * alpha * ik
    cmplx C_mfie;       // eta * (1-alpha)

    // File I/O
    int nsrcs;          
    double eleng;
    std::string lengStr;
};