#pragma once


#include <vector>


// Brine velocity coefficients
static const double VFI[] = {  1.166672   ,  1.448732   , -2.325958    };
static const double VGI[] = { -9.089008e-3,  3.746890e-5, -4.104040e-8 };
static const double VHI[] = { -2.453360e-4,  3.843340e-6, -2.108510e-5 };
static const double VII[] = { -6.475100e-4,  4.919670e-5,  0.          };
static const double AI[]  = {  2.118370e1 , -5.449600e2 , 2.760800e-1, -4.600000e-4, -5.801000e-2, 1.347200, -7.700000e-4 };
static const double BI[]  = { -4.555800   , -1.996500e2 , 5.495900e-1, -2.960000e-3, -9.281000e-2, 2.281300, -9.600000e-4 };

// Brine density coefficients
static const double D1I[] = { -7.92500e-5, -1.93000e-6, -3.42540e-4,  0.        ,  0.      };
static const double D2I[] = {  1.09980e-3, -2.87550e-3, -3.58190e-3, -7.28770e-1,  1.92016 };
static const double D3I[] = { -7.64020e-3,  3.69630e-2,  4.36083e-2, -3.33661e-1,  1.18569 };
static const double D4I[] = {  3.74600e-4, -3.32800e-4, -3.34600e-4,  0.        ,  0.      };
static const double E1I[] = {  0.        ,  0         ,  1.35300e-1,  0.        ,  0.      };
static const double F1I[] = { -1.40900   , -3.61000e-1, -2.53200e-1,  0.        ,  9.21600 };
static const double F2I[] = {           0,  5.61400   ,  4.67820   , -3.07000e-1,  2.60690 };
static const double F3I[] = { -1.12700e-1,  2.04700e-1, -4.52000e-2,  0.        ,  0.      };
static const double A2I[] = { -1.210000e-3, -1.125000e-1,  1.738500e-4, -9.360000e-9 };
static const double B2I[] = { -9.090000e-3,  9.859000e-2,  2.018980e-4,  1.665200e-7 };


namespace Brine {

    double velocity(double P, double T, double Sa, double NaCl, double KCl, double CaCl2);
    std::vector<double> velocity(const std::vector<double>& P, const std::vector<double>& T,
                                const std::vector<double>& Sa, const std::vector<double>& NaCl,
                                const std::vector<double>& KCl, const std::vector<double>& CaCl2);

    double density(double P, double T, double Sa, double NaCl, double KCl, double CaCl2);
    std::vector<double> density(const std::vector<double>& P, const std::vector<double>& T,
                               const std::vector<double>& Sa, const std::vector<double>& NaCl,
                               const std::vector<double>& KCl, const std::vector<double>& CaCl2);

    double bulk_modulus(double P, double T, double Sa, double NaCl, double KCl, double CaCl2);
    std::vector<double> bulk_modulus(const std::vector<double>& P, const std::vector<double>& T,
                                    const std::vector<double>& Sa, const std::vector<double>& NaCl,
                                    const std::vector<double>& KCl, const std::vector<double>& CaCl2);

    double viscosity(double T, double Sa);
    std::vector<double> viscosity(const std::vector<double>& T, const std::vector<double>& Sa);

    double resistivity(double T, double Sa);
    std::vector<double> resistivity(const std::vector<double>& T, const std::vector<double>& Sa);

    double gas_solubility(double P, double T, double Sa);
    std::vector<double> gas_solubility(const std::vector<double>& P, const std::vector<double>& T,
                                       const std::vector<double>& Sa);

}
