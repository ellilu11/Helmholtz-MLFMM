#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <type_traits>

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

template <typename T>
void importLines(const std::filesystem::path& path, std::vector<T>& vec) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Unable to find file");
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        T ele;

        if (iss >> ele) vec.emplace_back(ele);
        else throw std::runtime_error("Unable to parse line");
    }
}

void makeDir(const std::filesystem::path& dir) {
    std::error_code ec;
    if (std::filesystem::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";
}