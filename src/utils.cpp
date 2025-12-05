#include <functional>
#include <iostream>
#include <cmath>

#include "co2.hpp"
#include "utils.hpp"


// double co2_secant(double x0, double P, double T, double tol, double maxiter) {
//     if (tol == -1) tol = 1e-10;
//     if (maxiter == -1) maxiter = 100;
// 
//     std::function<double(double, double, double)> func = CO2::density_objective;
// 
//     double eps = 1e-4;
//     double x1 = x0 + eps;
//     double y0 = func(x0, P, T);
//     double y1 = func(x1, P, T);
// 
//     double x = x1 - y1*(x1-x0)/(y1-y0);
//     x0 = x1;
//     x1 = x;
// 
//     int iter = 1;
//     for (int i=0; i<maxiter; ++i) {
//         y0 = func(x0, P, T);
//         y1 = func(x1, P, T);
//         x = x1 - y1*(x1-x0)/(y1-y0);
//         x0 = x1;
//         x1 = x;
//         iter += 1;
//         if (std::abs(x1-x0) <= tol) break;
//     }
// 
//     if (iter > maxiter)
//         std::cerr << "Failed to converge." << std::endl;
// 
//     return x1;
// }


double co2_secant(double x0, double P, double T, double tol, int maxiter) {
    std::function<double(double,double,double)> func = CO2::density_objective;

    // Second initial guess
    double x1 = x0 * 1.0001;          // or better: user-supplied
    double f0 = func(x0, P, T);
    double f1 = func(x1, P, T);

    for (int i = 0; i < maxiter; ++i) {

        double denom = (f1 - f0);
        if (std::abs(denom) < 1e-14) {
            std::cerr << "Secant method: division by zero." << std::endl;
            return x1;
        }

        double x2 = x1 - f1 * (x1 - x0) / denom;

        // convergence check
        if (std::abs(x2 - x1) <= tol)
            return x2;

        // shift
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = func(x1, P, T);
    }

    std::cerr << "Secant method failed to converge." << std::endl;
    return x1;
}
