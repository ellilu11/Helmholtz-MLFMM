#include <filesystem>
#include <random>
#include "dipole.h"
#include "mesh/triangle.h"

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

SrcVec importSources(std::shared_ptr<Exc::PlaneWave> Einc)
{
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
    Mesh::Triangle::buildQuadCoeffs(config.quadPrec);

    const string configPath = "config/rwg/test/n"+to_string(config.nsrcs)+"adj/";
    // const string configPath = "config/rwg/sph"+to_string(config.nsrcs)+"/";
    auto srcs = Mesh::importMesh(
        configPath+"vertices.txt",
        configPath+"faces.txt",
        configPath+"rwgs.txt",
        Einc);

    // Mesh::refineMesh(srcs);

    return srcs;
}

