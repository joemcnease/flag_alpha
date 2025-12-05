#include <cmath>
#include <ctime>
#include <vector>
#include <stdexcept>

#include "license.hpp"
#include "utils.hpp"
#include "h2o.hpp"


namespace H2O {


    /*
     * Velocity of pure water.
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa)
     *   T : double
     *     Saturated vapor temperature (C).
     *
     * Returns:
     *   double
     *     Pure water velocity at corresponding (P, T).
     */
    double velocity(double P, double T) {
        double P2 = P*P;
        double P3 = P2*P;
        double av = AVI[0] + AVI[1]*P + AVI[2]*P2 + AVI[3]*P3;
        double bv = BVI[0] + BVI[1]*P + BVI[2]*P2 + BVI[3]*P3;
        double cv = CVI[0] + CVI[1]*P + CVI[2]*P2 + CVI[3]*P3;
        double dv = DVI[0] + DVI[1]*P + DVI[2]*P2 + DVI[3]*P3;
        double ev = EVI[0] + EVI[1]*P + EVI[2]*P2 + EVI[3]*P3;

        double T2 = T*T;
        double T3 = T2*T;
        double T4 = T3*T;

        double V = 0.001*(av + bv*T + cv*T2 + dv*T3 + ev*T4);

        return V;
    }


    /*
     * Velocity of pure water for an array of (P, T) values.
     *
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     */
    std::vector<double> velocity(const std::vector<double>& P, const std::vector<double>& T) {
        if (!are_same_size(P, T))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = velocity(P[i], T[i]);
        }

        return out;
    }


    /*
     * Density of pure water.
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Saturated vapor temperature (C).
     *
     * Returns:
     *   double
     *     Pure water density at corresponding (P, T).
     */
    double density(double P, double T) {
        double Tc = T/100.;
        double P0 = 70.;
        double Tc2 = Tc*Tc;
        double rhow0 = (DI[0]*Tc2 + DI[1]*Tc + DI[2])/(DI[3]*Tc2 + DI[4]*Tc + 1.);
        double Ew    = (EI[0]*Tc2 + EI[1]*Tc + EI[2])/(EI[3]*Tc2 + EI[4]*Tc + 1.);
        double Fw    = (FI[0]*Tc2 + FI[1]*Tc + FI[2])/(FI[3]*Tc2 + FI[4]*Tc + 1.);
        double Iw = std::log(std::abs(Ew*(P/P0) + Fw))/Ew;
        double Iw0 = std::log(std::abs(Ew + Fw))/Ew;

        double rhow = rhow0*std::exp(Iw - Iw0);

        return rhow;
    }


    /*
     * Density of pure water for an array of (P, T) values.
     *
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     */
    std::vector<double> density(const std::vector<double>& P, const std::vector<double>& T) {
        if (!are_same_size(P, T))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = density(P[i], T[i]);
        }

        return out;
    }


    /* Bulk modulus of pure water.
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Saturated vapor temperature (C).
     *
     * Returns:
     *   double
     *     Pure water bulk modulus at corresponding (P, T).
     */
    double bulk_modulus(double P, double T) {
        double rhow = density(P, T);
        double Vw = velocity(P, T);

        double Kw = rhow*Vw*Vw;

        return Kw;
    }


    /*
     * Bulk modulus of pure water for an array of (P, T) values.
     *
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     */
    std::vector<double> bulk_modulus(const std::vector<double>& P, const std::vector<double>& T) {
        if (!are_same_size(P, T))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = bulk_modulus(P[i], T[i]);
        }

        return out;
    }


    /* Viscosity of pure water.
     *
     * Parameters:
     *   T : double
     *     Temperature (C).
     *
     * Returns:
     *   double
     *     Pure water viscosity (centipose).
     */
    double viscosity(double T) {
        double eta = 0.1 + 1.65*std::exp(-0.057138*std::pow(T, 0.8));

        return eta;
    }


    /*
     * Viscosity of pure water for an array of T values.
     *
     * Parameters:
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     */
    std::vector<double> viscosity(const std::vector<double>& T) {
        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(T.size());
        for (size_t i=0; i<T.size(); ++i) {
            out[i] = viscosity(T[i]);
        }

        return out;
    }


    /* Saturated vapor pressure of pure water.
     *
     * Parameters:
     *   T : double
     *     Temperature (C).
     *
     * Returns:
     *   double
     *     Saturated vapor pressure (MPa).
     */
    double saturated_vapor_pressure(double T) {
        double A, B, C;
        if (T >= 0.01 && T < 100.) {
            A = 0.07131;
            B = 1730.;
            C = 233.426;
        } else if (T >= 100. && T <= 373.946) {
            A = 0.07131;
            B = 1730.;
            C = 233.426;
        }
        double exp = A-B/(C+T)-4.;
        double P = 1.33322 * std::pow(10., exp);

        return P;
    }


    /*
     * Saturated vapor pressure of pure water for an array of T values.
     *
     * Parameters:
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     */
    std::vector<double> saturated_vapor_pressure(const std::vector<double>& T) {
        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(T.size());
        for (size_t i=0; i<T.size(); ++i) {
            out[i] = saturated_vapor_pressure(T[i]);
        }

        return out;
    }


    /* Saturated vapor temperature of pure water.
     *
     * Parameters:
     *   P : double
     *     Saturated vapor pressure (MPa).
     *
     * Returns:
     *   double
     *     Saturated vapor temperature (C).
     */
    double saturated_vapor_temperature(double P) {
        double A, B, C;
        if (P >= 0.0006 && P < 0.1019) {
            A = 0.07131;
            B = 1730.;
            C = 233.426;
        } else if (P >= 0.1019 && P <= 21.7175) {
            A = 0.07131;
            B = 1730.;
            C = 233.426;
        }
        double T = B/(A-std::log10(7500.62*P)) - C;

        return T;
    }


    /*
     * Saturated vapor temperature of pure water for an array of P values.
     *
     * Parameters:
     *   T : double[]
     *     Array of saturated vapor pressure (MPa).
     */
    std::vector<double> saturated_vapor_temperature(const std::vector<double>& P) {
        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = saturated_vapor_temperature(P[i]);
        }

        return out;
    }


    /* Solubility of petroleum gas in pure water (GWR).
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Temperature (C).
     *
     *  Returns:
     *    double
     *      Solubility of petroleum gas in pure water at corresponding (P, T).
     */
    double gas_solubility(double P, double T) {
        double T2 = T*T;
        double T3 = T2*T;
        double T4 = T3*T;
        double A = 8.15839 - 6.12265e-2*T + 1.91663e-4*T2 - 2.1654e-7*T3;
        double B = 1.01021e-2 - 7.44241e-5*T + 3.05553e-7*T2 - 2.94883e-10*T3;
        double C = (-9.02505 + 0.130237*T - 8.53425e-4*T2 + 2.34122e-6*T3 - 2.37049e-9*T4)*1e-7;

        double Rsw = A + B*P + C*P*P;

        return Rsw;
    }


    /* Gas H2O for an array of (P, T) values.
     * 
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     */
    std::vector<double> gas_solubility(const std::vector<double>& P, const std::vector<double>& T) {
        if (!are_same_size(P, T))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = gas_solubility(P[i], T[i]);
        }

        return out;
    }


} // END NAMESPACE H2O


