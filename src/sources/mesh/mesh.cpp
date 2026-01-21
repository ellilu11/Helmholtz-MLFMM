#include "mesh.h"

void Mesh::importVertices(const std::filesystem::path& vpath) {
    std::ifstream file(vpath);
    if (!file) throw std::runtime_error("Unable to find file");
    std::string line;

    while (getline(file, line)) {
        std::istringstream iss(line);
        vec3d vertex;

        if (iss >> vertex) glVerts.push_back(vertex);
        else throw std::runtime_error("Unable to parse line");
    }

    nverts = glVerts.size(); // cache number of coarse vertices
}

void Mesh::importTriangles(const std::filesystem::path& tpath) {
    std::ifstream file(tpath);
    std::string line;
    if (!file) throw std::runtime_error("Unable to find file");

    int iTri = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        vec3i iVerts;

        if (iss >> iVerts) glTris.emplace_back(iVerts, iTri++);
        else throw std::runtime_error("Unable to parse line");
    }

    ntris = glTris.size(); // cache number of coarse triangles
}

SrcVec Mesh::importRWGs(
    const std::filesystem::path& rpath,
    const std::shared_ptr<Excitation::PlaneWave> Einc)
{
    std::ifstream file(rpath);
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
    const std::shared_ptr<Excitation::PlaneWave> Einc) 
{
    importVertices(vpath); 

    importTriangles(tpath);

    return importRWGs(rpath, std::move(Einc));
}

void Mesh::refineMesh(const SrcVec& rwgs) {
    Triangle::refineVertices();
    Triangle::buildSubtris();
    Triangle::buildEdgeToTri();

    SubRWG::buildSubRWGs();
    SubRWG::buildVertsToSubrwgs(nverts);

    for (const auto& rwg : rwgs)
        dynamic_pointer_cast<SrcRWG>(rwg)->buildSubIdx();
    //    dynamic_pointer_cast<SrcRWG>(rwg)->buildBC();
}