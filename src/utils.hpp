#pragma once


// double co2_secant(double x0, double P, double T, double tol, double maxiter);
double co2_secant(double x0, double P, double T, double tol = 1e-10, int maxiter = 100);

template <typename T0, typename... Ts>
bool are_same_size(T0 const& first, Ts const&... rest) {
    return ((first.size() == rest.size()) && ...);
}

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T>
T to_API(T rho) {
    return 141.5/rho - 131.5;
}

