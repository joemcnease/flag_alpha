#pragma once


#include <vector>


// Regression coefficients for CO2 equation of state. See Span and Wagner (1996).
static const double SWn[] = {  0.38856823203161   ,  0.29385475942740e1 , -0.55867188534934e1 ,
                              -0.76753199592477   ,  0.31729005580416   ,  0.54803315897767   ,
                               0.12279411220335   ,  0.21658961543220e1 ,  0.15841735109724e1 ,
                              -0.23132705405503   ,  0.58116916431436e-1, -0.55369137205382   ,
                               0.48946615909422   , -0.24275739843501e-1,  0.62494790501678e-1,
                              -0.12175860225246   , -0.37055685270086   , -0.16775879700426e-1,
                              -0.11960736637987   , -0.45619362508778e-1,  0.35612789270346e-1,
                              -0.74427727132052e-2, -0.17395704902432e-2, -0.21810121289527e-1,
                               0.24332166559236e-1, -0.37440133423463e-1,  0.14338715756878   ,
                              -0.13491969083286   , -0.23151225053480e-1,  0.12363125492901e-1,
                               0.21058321972940e-2, -0.33958519026368e-3,  0.55993651771592e-2,
                              -0.30335118055646e-3, -0.21365488688320e3 ,  0.26641569149272e5 ,
                              -0.24027212204557e5 , -0.28341603423999e3 ,  0.21247284400179e3 ,
                              -0.66642276540751   ,  0.72608632349897   ,  0.55068668612842e-1  };
static const double SWd[] = { 1, 1, 1, 1, 2, 2, 3, 1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7,
                             8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8, 2, 2, 2, 3, 3 };
static const double SWt[] = { 0.00, 0.75, 1.00, 2.00, 0.75, 2.00, 0.75, 1.50, 1.50, 2.50,
                             0.00, 1.50, 2.00, 0.00, 1.00, 2.00, 3.00, 6.00, 3.00, 6.00,
                             8.00, 6.00, 0.00, 7.00, 12.0, 16.0, 22.0, 24.0, 16.0, 24.0,
                             8.00, 2.00, 28.0, 14.0, 1.00, 0.00, 1.00, 3.00, 3.00 };
static const double SWc[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
                             2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6 };
static const double SWalpha[]  = { 25, 25, 25, 15, 20 };
static const double SWbeta[]   = { 325, 300, 300, 275, 275 };
static const double SWgamma[]  = { 1.16, 1.19, 1.19, 1.25, 1.22 };
static const double SWepsil[]  = { 1.00, 1.00, 1.00, 1.00, 1.00 };
static const double SWa[]      = { 3.5, 3.5, 3.0 };
static const double SWb[]      = { 0.875, 0.925, 0.875 };
static const double SWbeta2[]  = { 0.3, 0.3, 0.3 };
static const double SWA[]      = { 0.7, 0.7, 0.7 };
static const double SWB[]      = { 0.3, 0.3, 1.0 };
static const double SWC[]      = { 10.0, 10.0, 12.5 };
static const double SWD[]      = { 275, 275, 275 };
static const double SWao[]     = { 8.37304456, -3.70454304, 2.5, 1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678 };
static const double SWthetao[] = { 0, 0, 0, 3.15163, 6.11190, 6.77708, 11.32384, 27.08792 };


namespace CO2 {


    double density_objective(double x, double P, double T);

    double pressure(double rho, double T);
    std::vector<double> pressure(const std::vector<double>& rho, const std::vector<double>& T);

    double velocity(double P, double T);
    std::vector<double> velocity(const std::vector<double>& P, const std::vector<double>& T);

    double density(double P, double T);
    std::vector<double> density(const std::vector<double>& P, const std::vector<double>& T);

    double bulk_modulus(double P, double T);
    std::vector<double> bulk_modulus(const std::vector<double>& P, const std::vector<double>& T);

}
