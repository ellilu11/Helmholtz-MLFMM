#pragma once

#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <type_traits>

class Source;
namespace Mesh { class RWG; }
namespace Exct { class PlaneWave; }

void getDigit(std::istringstream& iss, char ch) {
    while (iss.get(ch)) {
        if (std::isdigit(static_cast<unsigned char>(ch))) {
            iss.unget();
            break;
        }
    }
}

template <typename Tin, typename Tout>
void importVec(const std::filesystem::path& path, std::vector<Tout>& vec) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Unable to find file");
    std::string line;

    size_t idx = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        Tin arg, arg1;
        double arg2;

        // Plane waves
        if constexpr (std::is_same_v<Tout, std::unique_ptr<Exct::PlaneWave>>) {
            if (iss >> arg >> arg1 >> arg2)
                vec.push_back(std::make_unique<Exct::PlaneWave>(arg, arg1, arg2));
            else throw std::runtime_error("Unable to parse line");
        } else if (iss >> arg) {
            // Vertices
            if constexpr (std::is_same_v<Tin, Tout>)
                vec.push_back(arg);
            // Sources
            else if constexpr (std::is_same_v<Tout, std::shared_ptr<Source>>)
                vec.push_back(std::make_shared<Mesh::RWG>(arg, idx++));
                // vec.push_back(std::make_shared<Dipole>(arg, arg1, idx++));
            // Triangles
            else
                vec.emplace_back(arg, idx++);
        } else throw std::runtime_error("Unable to parse line");
    }
}

void makeDir(const std::filesystem::path& dir) {
    std::error_code ec;
    if (std::filesystem::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";
}