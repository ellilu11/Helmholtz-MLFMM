#include "triangle.h"

using namespace std;

std::vector<vec3d> importVertices(const filesystem::path& fpath) {
    ifstream file(fpath);
    if (!file) throw std::runtime_error("Unable to find file");
    string line;
    std::vector<vec3d> vList;

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
            triangles.emplace_back(make_shared<Triangle>(vIdx,vList,prec));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return triangles;
}

int Triangle::quadPrec2Int(const Precision quadPrec) {
    return
        [&]() {
        switch (quadPrec) {
            case Precision::LOW:     return 1;
            case Precision::MEDIUM:  return 3;
            case Precision::HIGH:    return 7;
        };
        } ();
}

Triangle::Triangle(
    const vec3i& vIdx,
    const std::vector<vec3d>& vList,
    const Precision quadPrec)
    : vIdx(vIdx),
    vertices({ vList[vIdx[0]], vList[vIdx[1]], vList[vIdx[2]] }),
    center((vertices[0] + vertices[1] + vertices[2]) / 3.0)
{
    buildQuads(quadPrec);
};

void Triangle::buildQuads(const Precision quadPrec) {

    switch (quadPrec) {
        case Precision::LOW :
            quadNodes.push_back(center);
            quadWeight = 1.0/2.0;
            break;

        case Precision::MEDIUM :
            quadNodes.resize(3);
            quadNodes[0] = 2.0/3.0*vertices[0] + 1.0/6.0*vertices[1] + 1.0/6.0*vertices[2];
            quadNodes[1] = 1.0/6.0*vertices[0] + 2.0/3.0*vertices[1] + 1.0/6.0*vertices[2];
            quadNodes[2] = 1.0/6.0*vertices[0] + 1.0/6.0*vertices[1] + 2.0/3.0*vertices[2];
            quadWeight = 1.0/6.0;
            break;

        case Precision::HIGH :
            quadNodes.resize(7);
            // TODO: 7-point quadrature
            break;
    }
}