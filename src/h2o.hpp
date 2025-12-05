#pragma once


#include <vector>


// Pure water velocity coefficients
static const double AVI[] = {  1.408300e3 ,  1.593621   ,  1.707648e-3 , -4.392520e-6  };
static const double BVI[] = {  4.440552   , -2.830638e-3, -3.115400e-5 ,  5.247150e-8  };
static const double CVI[] = { -4.072162e-2,  1.489160e-4, -4.361920e-7 ,  8.588850e-10 };
static const double DVI[] = {  1.145250e-4, -8.154600e-7,  3.996740e-9 , -8.201520e-12 };
static const double EVI[] = { -1.579890e-7,  1.801380e-9, -1.019810e-11,  2.056540e-14 };

// Pure water density coefficients
static const double DI[] = { -1.272130e-1,  6.454860e-1,  1.032650  , -7.029100e-2,  6.395890e-1 };
static const double EI[] = {  4.221000   , -3.478000   ,  6.221000  ,  5.182000e-1, -4.405000e-1 };
static const double FI[] = { -1.140300e1 ,  2.993200e1 ,  2.795200e1,  2.068400e-1,  3.768000e-1 };


// Fluid
namespace H2O {


    double velocity(double P, double T);
    std::vector<double> velocity(const std::vector<double>& P, const std::vector<double>& T);

    double density(double P, double T);
    std::vector<double> density(const std::vector<double>& P, const std::vector<double>& T);

    double bulk_modulus(double P, double T);
    std::vector<double> bulk_modulus(const std::vector<double>& P, const std::vector<double>& T);

    double viscosity(double T);
    std::vector<double> viscosity(const std::vector<double>& T);

    double saturated_vapor_pressure(double T);
    std::vector<double> saturated_vapor_pressure(const std::vector<double>& T);

    double saturated_vapor_temperature(double P);
    std::vector<double> saturated_vapor_temperature(const std::vector<double>& P);

    double gas_solubility(double P, double T);
    std::vector<double> gas_solubility(const std::vector<double>& P, const std::vector<double>& T);

} // END NAMESPACE H2O
