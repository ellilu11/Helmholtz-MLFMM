#include <filesystem>
#include <random>
#include "sources/dipole.h"
#include "sources/mesh/mesh.h"

using namespace std; // TODO: Remove

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

pair<SrcVec, shared_ptr<Excitation::PlaneWave>> importFromConfig(const Config& config) 
{
    auto Einc = Excitation::importPlaneWave("config/pwave.txt");

    /* Dipole sources
    const auto fpath = makePath(config);
    SrcVec srcs;
    switch (config.mode) {
        case Mode::READ:
            srcs = importDipoles(fpath, Einc);
            break;
            
        case Mode::WRITE: {
            srcs = makeDipoles<uniform_real_distribution<double>>(config, Einc);

            ofstream srcFile(fpath);
            for (const auto& src : srcs) srcFile << *(dynamic_pointer_cast<Dipole>(src));
            break;
        }
    }
    cout << "   Source file:     " << fpath.generic_string() << '\n';
    */

    // RWG sources
    Mesh::Triangle::buildQuadCoeffs(config.quadPrec); // build triangle quadrature coeffs

    const string configPath = "config/rwg/test/n"+to_string(config.nsrcs)+"adj/";
    // const string configPath = "config/rwg/sph"+to_string(config.nsrcs)+"/";
    auto srcs = Mesh::importMesh(
        configPath+"vertices.txt",
        configPath+"faces.txt",
        configPath+"rwgs.txt",
        Einc);
    // Mesh::refineMesh(srcs);
    //

    cout << fixed << setprecision(3);
    cout << "   Mode:            " << (config.mode == Mode::READ ? "READ" : "WRITE") << '\n';
    cout << "   # Sources:       " << srcs.size() << '\n';
    cout << "   Digit precision: " << config.digits << '\n';
    cout << "   Interp order:    " << config.interpOrder << '\n';
    cout << "   RWG quad rule:   " << Mesh::Triangle::getNumQuads() << "-point\n";
    cout << "   Max node RWGs:   " << config.maxNodeSrcs << '\n';
    cout << "   Root length:     " << config.rootLeng << '\n';
    cout << "   Wave number:     " << Einc->wavenum << "\n\n";

    return make_pair(srcs, Einc);
}

