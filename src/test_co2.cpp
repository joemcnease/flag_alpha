#include <iostream>
#include <vector>

#include "co2.hpp"


int main(int argc, char *argv[]) {
    double P = 40;
    double T = 100;
    double V = CO2::velocity(P, T);
    std::cout << V << std::endl;

    return 0;
}
