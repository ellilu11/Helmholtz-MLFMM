#include "mesh.h"
#include "srcrwg.h"

void Mesh::importVertices(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Unable to find file");
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        vec3d vertex;

        if (iss >> vertex) glVerts.push_back(vertex);
        else throw std::runtime_error("Unable to parse line");
    }

    nverts = glVerts.size(); // record number of coarse vertices
}

void Mesh::importTriangles(const std::filesystem::path& path) {
    std::ifstream file(path);
    std::string line;
    if (!file) throw std::runtime_error("Unable to find file");

    int iTri = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        vec3i iVerts;

        if (iss >> iVerts) glTris.emplace_back(iVerts, iTri++);
        else throw std::runtime_error("Unable to parse line");
    }

    ntris = glTris.size(); // record number of coarse triangles
    triToRWGs.resize(ntris);
}

SrcVec Mesh::importRWGs(
    const std::filesystem::path& path,
    const std::shared_ptr<Exc::PlaneWave> Einc)
{
    std::ifstream file(path);
    std::string line;
    if (!file) throw std::runtime_error("Unable to find file");
    SrcVec rwgs;
    size_t iSrc = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        Eigen::Vector4i idx4;

        if (iss >> idx4)
            rwgs.push_back(std::make_shared<SrcRWG>(Einc, iSrc++, idx4));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return rwgs;
}

SrcVec Mesh::importMesh(
    const std::filesystem::path& vpath,
    const std::filesystem::path& tpath,
    const std::filesystem::path& rpath,
    std::shared_ptr<Exc::PlaneWave> Einc) 
{
    importVertices(vpath); 

    importTriangles(tpath);

    return importRWGs(rpath, std::move(Einc));
}

void Mesh::printScattered(const SrcVec& srcs, const std::string& fname, int nth) {
    namespace fs = std::filesystem;
    fs::path dir = "out/ff/pz_kx_r5.0";
    std::error_code ec;

    std::cout << " Computing scattered farfield...\n";

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    std::ofstream farfile(dir/fname);
    farfile << std::setprecision(15) << std::scientific;

    // Also print out angles (coordinates of farsols)
    std::ofstream thfile(dir/"thetas.txt"); // phfile(dir/"phis.txt");
    thfile << std::setprecision(15) << std::scientific;

    double phi = 0.0; // pick phi = 0
    double rcsSum = 0.0;
    for (int ith = 0; ith < nth; ++ith) {
        double theta = (ith+0.5)*PI/static_cast<double>(nth);
        const vec3d& rhat = Math::fromSph(vec3d(1.0, theta, phi));

        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += states.currents[src->getIdx()] * src->getFarAlongDir(k*rhat);

        const vec2cd& far = Phys::C * k * Math::toThPh(theta, phi) * dirFar;

        farfile << norm(far[0]) << ' ' << norm(far[1]) << '\n'; // squared magnitude of theta and phi components
        thfile << theta << '\n';

        rcsSum += norm(far[0]) + norm(far[1]);
    }

    std::cout << " Average RCS: " << rcsSum/nth << "\n";
}

/*
void Mesh::printRWGs(const SrcVec& rwgs, const std::string& fname) {
    std::ofstream rwgfile(dir/fname);

    for (const auto& rwg : rwgs) {
        const auto& vertsC = dynamic_pointer_cast<RWG>(rwg)->getVertsC();
        const vec3d center = (vertsC[0] + vertsC[1]) / 2.0;

        rwgfile << center << '\n';
    }
}*/