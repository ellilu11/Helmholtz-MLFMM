#pragma once

#include <array>
#include <complex>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

using cmplx = std::complex<double>;

// TODO: unalias these?
using pair2i = std::pair<int, int>;
using pair2d = std::pair<double, double>;
using pair2di = std::pair<double, int>;
using pair2cd = std::pair<cmplx, cmplx>;

std::ostream& operator<< (std::ostream& os, cmplx z) {
    //char sign = z.imag() >= 0.0 ? '+' : '-';
    //os << z.real() << sign << abs(z.imag()) << 'i';
    os << z.real() << ' ' << z.imag();
    return os;
}

template <typename T>
std::vector<T> operator+= (std::vector<T>& lhs, const std::vector<T>& rhs) {
    assert(lhs.size() == rhs.size());
    for (size_t i = 0; i < lhs.size(); ++i)
        lhs[i] += rhs[i];
    return lhs;
}

template <typename T>
std::vector<T> operator+ (std::vector<T> lhs, const std::vector<T>& rhs) {
    lhs += rhs;
    return lhs;
}

template <typename T>
std::vector<T> operator*= (std::vector<T>& vec, T val) {
    for (size_t i = 0; i < vec.size(); ++i)
        vec[i] *= val;
    return vec;
}

template <typename T>
std::vector<T> operator* (T val, std::vector<T> vec) {
    vec *= val;
    return vec;
}

template <typename T>
std::ostream& operator<< (std::ostream& os, const std::vector<T>& vec) {
    for (const auto& val : vec) os << val << " ";
    os << '\n';
    return os;
}

