#include "rwg.h"

using namespace std;

RWGVec importRWG(
    const filesystem::path& fpath, 
    const std::vector<vec3d>& vertices,
    const TriVec& triangles,
    const shared_ptr<Src> Einc)
{
    ifstream file(fpath);
    string line;
    if (!file) throw std::runtime_error("Unable to find file");
    RWGVec rwgs;

    while (getline(file, line)) {
        istringstream iss(line);
        Eigen::Vector4i idx;

        if (iss >> idx)
            rwgs.emplace_back(make_shared<RWG>(idx, vertices, triangles, Einc));
        else
            throw std::runtime_error("Unable to parse line");
    }

    return rwgs;
}

RWG::RWG(const Eigen::Vector4i& idx,
    const std::vector<vec3d>& vertices,
    const TriVec& triangles,
    const std::shared_ptr<Src> Einc)
    : v0(vertices[idx[0]]), v1(vertices[idx[1]]),
      triPlus(triangles[idx[2]]), triMinus(triangles[idx[3]]),
      vCenter((v0+v1)/2.0),
      Einc(Einc)
{
    for (const auto& vIdx : triPlus->vIdx)
        if (vIdx != idx[0] && vIdx != idx[1])
            vPlus = vertices[vIdx];

    for (const auto& vIdx : triMinus->vIdx)
        if (vIdx != idx[0] && vIdx != idx[1])
            vMinus = vertices[vIdx];

    leng = (v0-v1).norm();

    buildCurrent();

    // std::cout << '(' << v0 << ") (" << v1 << ") ("
    //    << vPlus << ") (" << vMinus << ") " << leng << '\n';
};

void RWG::buildRHS() {

    cmplx rhs(0,0);

    // TODO: Rewrite this
    /*std::vector<vec3d> midPlus = { vCenter, (v0 + vPlus)/2.0, (v1 + vPlus)/2.0 };
    std::vector<vec3d> midMinus = { vCenter, (v0 + vMinus)/2.0, (v1 + vMinus)/2.0 };

    for (const auto& node : midPlus)
       rhs += -leng / (6.0 * triPlus->getArea()) * 
            (node - vPlus).dot(Einc->amplitude * Einc->pol) * 
            Math::expI(Einc->wavevec.dot(node));

    for (const auto& node : midMinus)
        rhs -= -leng / (6.0 * triPlus->getArea()) *
            (node - vPlus).dot(Einc->amplitude * Einc->pol) *
            Math::expI(Einc->wavevec.dot(node));
    */
    
}

void RWG::buildCurrent() {
    current = 1.0;
}
