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

void Mesh::getScattered(
    const SrcVec& srcs, const std::filesystem::path& dir, const std::string& fname, int nth) 
{
    makeDir(dir);
    std::ofstream farfile(dir/fname), thfile(dir/"angles.txt");
    farfile << std::setprecision(15) << std::scientific;

    std::cout << " Computing scattered farfield...\n";

    double k = config.k;
    double rcsSum = 0.0;
    for (int ith = 0; ith < nth; ++ith) { // 2*nth to cover great circle
        double theta0 = (ith+0.5)*PI/static_cast<double>(nth); // in [0, 2*pi]
        double theta = (theta0 < PI) ? theta0 : 2*PI - theta0; // fold back to [0, pi]
        double phi = (theta0 < PI) ? 0.0 : PI; // fold back to [0, 2pi]
        assert(theta >= 0.0 && theta <= PI);

        vec3d rhat = Math::fromSph(vec3d(1.0, theta, phi));

        vec3cd dirFar = vec3cd::Zero();
        for (const auto& src : srcs)
            dirFar += Solver::currents[src->getIdx()] * src->getFarAlongDir(k*rhat);

        const vec2cd& far = Phys::C * k * Math::toThPh(theta, phi) * dirFar;

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