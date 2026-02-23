#include <iostream>
#include "../src/types.h"
// #include "../src/fmm/coeffs.h"

int main() {
    std::cout << "Testing std::vector arithmetic operators...\n";

    double a = 2.0;
    std::vector<double> v = { 1.0, 2.0, 3.0 };
    
    std::vector<double> w = a*v;
    std::cout << "v = " << v[0] << " " << v[1] << " " << v[2] << "\n";
    std::cout << "a*v = " << w[0] << " " << w[1] << " " << w[2] << "\n";

    v *= a;
    std::cout << "a*v = " << v[0] << " " << v[1] << " " << v[2] << "\n";

    /*
    std::cout << "Testing Coeffs arithmetic operators...\n";
    double a = 2.0;
    FMM::Coeffs C(3, 1.0);
    FMM::Coeffs D = a * C;
    std::cout << "C =\n" << C << "\n";
    std::cout << "a*C =\n" << a*C << "\n";
    std::cout << "a*C =\n" << D << "\n";
    std::cout << "C+D =\n" << C+D << "\n";
    */

    return 0;
}