#include "rwg.h"
#include "triangle.h"

SrcVec Mesh::importMesh(const std::filesystem::path& path)
{
    std::cout << " Importing mesh...\n";

    Triangle::buildQuadCoeffs(config.quadPrec);

    SrcVec rwgs;
    importVec<vec3d>(path/"vertices.txt", glVerts);
    importVec<vec3i>(path/"faces.txt", glTris);
    importVec<vec4i>(path/"rwgs.txt", rwgs);

    buildRootCoords();
    printNormals("out/nhats.txt");

    std::cout << "   # Sources:       " << rwgs.size() << '\n';
    std::cout << "   Root length:     " << rootLeng << " m\n";
    std::cout << "   Root center:     " << rootCenter << "\n\n";

    return rwgs;
}

void Mesh::buildRootCoords() {
    vec3d maxVert = vec3d::Constant(-std::numeric_limits<double>::infinity());
    vec3d minVert = vec3d::Constant(std::numeric_limits<double>::infinity());

    for (const auto& vert : glVerts) {
        maxVert = maxVert.cwiseMax(vert);
        minVert = minVert.cwiseMin(vert);
    }

    rootLeng = (maxVert - minVert).lpNorm<Eigen::Infinity>()
        * (1.0 + 1e-3); // add wiggle room
    rootCenter = 0.5*(maxVert + minVert);
}

void Mesh::printNormals(const std::string& fname) {
    std::ofstream file(fname);
    file << std::setprecision(15) << std::scientific;
    for (const auto& tri : glTris) {
        vec3d center = tri.getCenter();
        vec3d nhat = tri.getNormal();
        file << center[0] << ' ' << center[1] << ' ' << center[2] << ' ' 
            << nhat[0] << ' ' << nhat[1] << ' ' << nhat[2] << ' '
            << center.dot(nhat)
            << '\n';
    }
}