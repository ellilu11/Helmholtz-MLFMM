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

template <typename T>
std::ifstream& operator>>(std::ifstream& is, T& val) {
    std::string line;
    if (std::getline(is, line)) {
        std::istringstream iss(line);

        char ch = '\0';
        getDigit(iss, ch);

        if constexpr (std::is_enum_v<T>) {
            typename std::underlying_type<T>::type eval;

            while (iss >> eval) {
                val = static_cast<T>(eval);
                getDigit(iss, ch);
            }

        } else
            while (iss >> val)
                getDigit(iss, ch);
    }

    return is;
}

template <typename Tin, typename Tout>
void importVec(const std::filesystem::path& path, std::vector<Tout>& vec) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Unable to find file");
    std::string line;

    size_t idx = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        Tin ele, ele1;
        double ele2;

        // PlaneWave
        if constexpr (std::is_same_v<Tout, std::unique_ptr<Exct::PlaneWave>>) {
            if (iss >> ele >> ele1 >> ele2)
                vec.push_back(std::make_unique<Exct::PlaneWave>(ele, ele1, ele2));
            else throw std::runtime_error("Unable to parse line");
        } else if (iss >> ele) {
            // Vertices
            if constexpr (std::is_same_v<Tin, Tout>)
                vec.push_back(ele);
            // Sources
            else if constexpr (std::is_same_v<Tout, std::shared_ptr<Source>>)
                vec.push_back(std::make_shared<Mesh::RWG>(ele, idx++));
            // TODO: Handle Exct::PlaneWave and other source types
        // Triangles
            else
                vec.emplace_back(ele, idx++);
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