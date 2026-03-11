#include "rwg.h"
#include "triangle.h"

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

    // std::cout << " Sphere radius: " << glVerts[0].norm() << '\n';

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
    const std::filesystem::path& path, std::shared_ptr<Exc::PlaneWave> Einc)
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
            rwgs.push_back(std::make_shared<RWG>(Einc, iSrc++, idx4));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return rwgs;
}

void Mesh::buildRootCoords() {
    vec3d maxVert = vec3d::Constant(-std::numeric_limits<double>::infinity());
    vec3d minVert = vec3d::Constant(std::numeric_limits<double>::infinity());

    for (const auto& vert : glVerts) {
        maxVert = maxVert.cwiseMax(vert);
        minVert = minVert.cwiseMin(vert);
    }

    rootCenter = 0.5*(maxVert + minVert);
    rootLeng = (maxVert - minVert).lpNorm<Eigen::Infinity>() 
        * (1.0 + 1e-3); // add wiggle room
    std::cout << "   Root center: " << rootCenter << '\n';
    std::cout << "   Root length: " << rootLeng << " m\n\n";
}

SrcVec Mesh::importMesh(
    const std::filesystem::path& path, std::shared_ptr<Exc::PlaneWave> Einc) 
{
    std::cout << " Importing mesh...\n";

    Triangle::buildQuadCoeffs(config.quadPrec);

    importVertices(path/"vertices.txt");

    buildRootCoords();

    importTriangles(path/"faces.txt");

    return importRWGs(path/"rwgs.txt", std::move(Einc));
}

void Mesh::printScattered(
    const SrcVec& srcs, const std::filesystem::path& dir, const std::string& fname, int nth) {
    namespace fs = std::filesystem;
    std::error_code ec;

    std::cout << " Computing scattered farfield...\n";

    if (fs::create_directory(dir, ec))
        std::cout << "   Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << "   Error creating directory " << ec.message() << "\n";

    std::ofstream farfile(dir/fname);
    farfile << std::setprecision(15) << std::scientific;

    // Also print out angles (coordinates of farsols)
    std::ofstream thfile(dir/"angles.txt");

    double k = config.k;
    double rcsSum = 0.0;
    for (int ith = 0; ith < 2*nth; ++ith) {
        double theta0 = (ith+0.5)*PI/static_cast<double>(nth); // in [0, 2*pi]
        double theta = (theta0 < PI) ? theta0 : 2*PI - theta0; // fold back to [0, pi]
        double phi = (theta0 < PI) ? 0.0 : PI; // fold back to [0, 2pi]
        assert(theta >= 0.0 && theta <= PI);

        const vec3d& rhat = Math::fromSph(vec3d(1.0, theta, phi));

        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += Solver::currents[src->getIdx()] * src->getFarAlongDir(k*rhat);

        const vec2cd& far = iu*k*Phys::eta * Math::toThPh(theta, phi) * dirFar;

        double rcs = 4.0*PI/(k*k) * far.squaredNorm();
        rcsSum += rcs;
        farfile << rcs << '\n'; // squared magnitude of theta and phi components

        thfile << theta0 << ' ' << theta << ' ' << phi << '\n';
    }

    std::cout << "   Mean RCS: " 
        << std::setprecision(9) << rcsSum/nth << std::setprecision(3) << "\n";
    std::cout << "   File: " << (dir/fname).generic_string() << "\n\n";
}

void Mesh::printNormals(const std::string& fname) {
    std::ofstream file(fname);
    file << std::setprecision(15) << std::scientific;
    for (const auto& tri : glTris)
        file << tri.getCenter() << ' ' << tri.getNormal()
            << ' ' << tri.getCenter().dot(tri.getNormal())
            << '\n';
}