#include "triangle.h"

using namespace std;

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

int Triangle::prec2Int(const Precision quadPrec) {
    return
        [&]() {
        switch (quadPrec) {
            case Precision::LOW:     return 1;
            case Precision::MEDIUM:  return 3;
            case Precision::HIGH:    return 7;
        };
        } ();
}

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