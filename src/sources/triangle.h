#pragma once

class Triangle;

using TriVec = std::vector<std::shared_ptr<Triangle>>;

class Triangle {
public:
    Triangle()
        : vIdx(vec3i::Zero()),
          vertices( {zeroVec,zeroVec,zeroVec}), 
          center(zeroVec), area(0.0) { };

    Triangle(
        const vec3i&,
        const std::vector<vec3d>&,
        const Precision);

    vec3i getVidx() { return vIdx; }

    std::array<vec3d, 3> getVertices() { return vertices; }

    vec3d getCenter() const { return center; }

    double getArea() const { return area; }

    std::pair<std::vector<vec3d>, double> getQuads() {
        return std::make_pair(quadNodes, quadWeight);
    }

    static int prec2Int(const Precision);

    void buildQuads(const Precision);

private:
    vec3i vIdx;
    std::array<vec3d,3> vertices;
    vec3d center;
    double area;

    std::vector<vec3d> quadNodes;
    double quadWeight;
};