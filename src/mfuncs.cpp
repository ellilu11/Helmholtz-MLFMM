#include "MLFMA.h"

using namespace std;

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

vector<vec3d> importVertices(const filesystem::path& fpath) {
    ifstream file(fpath);
    if (!file) throw std::runtime_error("Unable to find file");
    string line;
    vector<vec3d> vList;

    while (getline(file, line)) {
        istringstream iss(line);
        vec3d vertex;

        if (iss >> vertex)
            vList.push_back(vertex);
        else
            throw std::runtime_error("Unable to parse line");
    }

    return vList;
}

TriVec importTriangles(
    const filesystem::path& fpath, const std::vector<vec3d>& vList, const Precision prec)
{
    ifstream file(fpath);
    string line;
    if (!file) throw std::runtime_error("Unable to find file");
    TriVec triangles;

    while (getline(file, line)) {
        istringstream iss(line);
        vec3i vIdx;

        if (iss >> vIdx)
            triangles.emplace_back(make_shared<Triangle>(vIdx, vList, prec));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return triangles;
}

SrcVec importRWG(
    const filesystem::path& vpath,
    const filesystem::path& tpath,
    const filesystem::path& rpath,
    const Precision quadPrec,
    const shared_ptr<PlaneWave> Einc)
{
    auto vertices = importVertices(vpath);

    auto triangles = importTriangles(tpath, vertices, quadPrec);

    ifstream file(rpath);
    string line;
    if (!file) throw std::runtime_error("Unable to find file");
    SrcVec rwgs;

    while (getline(file, line)) {
        istringstream iss(line);
        Eigen::Vector4i idx;

        if (iss >> idx)
            rwgs.emplace_back(make_shared<RWG>(Einc, idx, vertices, triangles));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return rwgs;
}

SrcVec importConfig(const Config& config, const string& configPath) {
    cout << " Importing config...\n";

    const auto fpath = makePath(config);

    shared_ptr<PlaneWave> Einc = make_shared<PlaneWave>();

    SrcVec srcs;

    switch (config.mode) {
        case Mode::READ:
            srcs = importDipoles(fpath);
            break;

        case Mode::WRITE: {
            srcs = makeDipoles<uniform_real_distribution<double>>(config, Einc);

            ofstream srcFile(fpath);
            for (const auto& src : srcs) srcFile << *src;
            break;
        }

    //const string configPath = "config/n480/";
    //auto srcs = importRWG(configPath+"vertices.txt",
    //                      configPath+"faces.txt",
    //                      configPath+"rwgs.txt",
    //                      config.quadPrec,
    //                      Einc);

    cout << "   # Sources:       " << srcs.size() << '\n';
    cout << "   RWG quad rule:   " << Triangle::prec2Int(config.quadPrec) << "-point\n";
    cout << "   Digit precision: " << config.digits << '\n';
    cout << "   Interp order:    " << config.interpOrder << '\n';
    cout << "   Max node RWGs:   " << config.maxNodeSrcs << '\n';
    cout << fixed << setprecision(3);
    cout << "   Root length:     " << config.rootLeng << '\n';
    cout << "   Wave number:     " << Einc->wavenum << "\n\n";

    return srcs;
}
