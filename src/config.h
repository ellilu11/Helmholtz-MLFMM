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
    Config(const std::string& fileName) {
        std::ifstream is(fileName);
        is >> mode >> quadPrec >> alpha >> k >> nsrcs 
           >> maxNodeSrcs >> digits >> interpOrder >> overInterp
           >> epsIter >> maxIter >> iluTol >> iluFactor
           >> eleng;

        if (alpha < 0.0 || alpha > 1.0)
            throw std::runtime_error("Alpha must be in [0,1]");
        ie = (alpha == 0.0) ? IE::MFIE : (alpha == 1.0) ? IE::EFIE : IE::CFIE;
        C_efie = -Phys::eta * alpha * iu * k;
        C_mfie = Phys::eta * (1.0 - alpha);

        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << eleng;
        lengStr = ss.str();
        
        std::cout << " *************************** \n";
        std::cout << " ***** Helmholtz-MLFMM ***** \n";
        std::cout << " *************************** \n";

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "   Mode:            " << getModeStr(mode) << '\n';
        std::cout << "   IE (alpha):      " << getIEStr(ie) << " (" << alpha << ")\n";
        std::cout << "   Wavenumber:      " << k << " /m\n";
        std::cout << "   Max in node:     " << maxNodeSrcs << '\n';
        std::cout << "   Digit precision: " << digits << '\n';
        std::cout << "   Interp order:    " << interpOrder << '\n';
        std::cout << "   Overinterp:      " << overInterp << '\n';
        std::cout << "   Tri quad rule:   " << getNumQuads(quadPrec) << "-point\n";
        std::cout << "   ILU drop tol:    " << iluTol << '\n';
        std::cout << "   ILU fill factor: " << iluFactor << "\n\n";
    }

    // General
    Mode mode;          // FMM, direct, or FMM+direct
    Precision quadPrec; // Triangle quadrature precision
    double alpha;       // CFIE weight
    double k;           // Wavenumber
    // double leps;     // TODO: Length error tolerance for nearfield
    int nsrcs;          // # sources (for file I/O only)

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

    // Elctrical length (for file I/O)
    double eleng;
    std::string lengStr;
};