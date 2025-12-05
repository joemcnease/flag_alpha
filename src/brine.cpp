#include <cmath>
#include <ctime>
#include <vector>
#include <stdexcept>

#include "license.hpp"
#include "utils.hpp"
#include "h2o.hpp"
#include "brine.hpp"


namespace Brine {


    /* Velocity of brine.
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Saturated vapor temperature (C).
     *   Sa : double
     *     Salinity (ppm).
     *   NaCl : double
     *     Weight percentage of NaCl of total salt in brine.
     *   KCl : double
     *     Weight percentage of KCl of total salt in brine.
     *   CaCl2 : double
     *     Weight percentage of CaCl2 of total salt in brine.
     *
     * Returns:
     *   double
     *     Brine velocity at corresponding (P, T, Sa, NaCl, KCl, CaCl2).
     */
    double velocity(double P, double T, double Sa, double NaCl, double KCl, double CaCl2) {
        double salt = Sa*1e-6;
        double salt1p5 = std::pow(salt, 1.5);
        double salt2 = salt*salt;
        double P2 = P*P;
        double T2 = T*T;
        double T3 = T2*T;

        double Vbs = VFI[0]*salt + VFI[1]*salt1p5 + VFI[2]*salt2;
        double Vbt = (VGI[0]*T + VGI[1]*T2 + VGI[2]*T3)*salt;
        double Vbp1 = (VHI[0]*P + VHI[1]*T*P + VHI[2]*P2)*salt;
        double Vbp2 = (VII[0]*P + VII[1]*P2)*salt1p5;
        double Vba = Vbs + Vbt + Vbp1 + Vbp2;
        double Vw = H2O::velocity(P, T);

        double Vbna = Vw + Vba;
        double Vbk = Vbna + 0.001*(AI[0] + AI[1]*salt + AI[2]*T + AI[3]*T2 + AI[4]*P + AI[5]*salt*T + AI[6]*salt*P);
        double Vbca = Vbna + 0.001*(BI[0] + BI[1]*salt + BI[2]*T + BI[3]*T2 + BI[4]*P + BI[5]*salt*T + BI[6]*salt*P);

        double Vb = 0.01*(NaCl*Vbna + KCl*Vbk + CaCl2*Vbca);

        return Vb;
    }


    /* Velocity of brine for an array of (P, T, Sa, NaCl, KCl, CaCl2) values.
     * 
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     *   Sa : double[]
     *     Array of salinity (ppm).
     *   NaCl : double[]
     *     Array of weight percentage of NaCl of total salt in brine.
     *   KCl : double[]
     *     Array of weight percentage of KCl of total salt in brine.
     *   CaCl2 : double[]
     *     Array of weight percentage of CaCl2 of total salt in brine.
     */
    std::vector<double> velocity(const std::vector<double>& P, const std::vector<double>& T,
                                const std::vector<double>& Sa, const std::vector<double>& NaCl,
                                const std::vector<double>& KCl, const std::vector<double>& CaCl2) {
        if (!are_same_size(P, T, Sa, NaCl, KCl, CaCl2))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = velocity(P[i], T[i], Sa[i], NaCl[i], KCl[i], CaCl2[i]);
        }

        return out;
    }


    /* Density of brine.
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Saturated vapor temperature (C).
     *   Sa : double
     *     Salinity (ppm).
     *   NaCl : double
     *     Weight percentage of NaCl of total salt in brine.
     *   KCl : double
     *     Weight percentage of KCl of total salt in brine.
     *   CaCl2 : double
     *     Weight percentage of CaCl2 of total salt in brine.
     *
     * Returns:
     *   double
     *     Brine density at corresponding (P, T, Sa, NaCl, KCl, CaCl2).
     */
    double density(double P, double T, double Sa, double NaCl, double KCl, double CaCl2) {
        double snacl = Sa/58.44 * 1e3/(1e6-Sa);
        double snaclsr = std::sqrt(snacl);
        double snacl1p5 = std::pow(snacl, 1.5);
        double snacl2 = snacl*snacl;

        double salt = Sa*1e-6;
        double T2 = T*T;
        double Tc = T/100;
        double Tc2 = Tc*Tc;
        double P0 = 70;

        double D1 = (D1I[0]*Tc2 + D1I[1]*Tc + D1I[2])/(D1I[3]*Tc2 + D1I[4]*Tc + 1.);
        double D2 = (D2I[0]*Tc2 + D2I[1]*Tc + D2I[2])/(D2I[3]*Tc2 + D2I[4]*Tc + 1.);
        double D3 = (D3I[0]*Tc2 + D3I[1]*Tc + D3I[2])/(D3I[3]*Tc2 + D3I[4]*Tc + 1.);
        double D4 = (D4I[0]*Tc2 + D4I[1]*Tc + D4I[2])/(D4I[3]*Tc2 + D4I[4]*Tc + 1.);
        double E1 = (E1I[0]*Tc2 + E1I[1]*Tc + E1I[2])/(E1I[3]*Tc2 + E1I[4]*Tc + 1.);
        double F1 = (F1I[0]*Tc2 + F1I[1]*Tc + F1I[2])/(F1I[3]*Tc2 + F1I[4]*Tc + 1.);
        double F2 = (F2I[0]*Tc2 + F2I[1]*Tc + F2I[2])/(F2I[3]*Tc2 + F2I[4]*Tc + 1.);
        double F3 = (F3I[0]*Tc2 + F3I[1]*Tc + F3I[2])/(F3I[3]*Tc2 + F3I[4]*Tc + 1.);
        double Ew = (EI[0]*Tc2 + EI[1]*Tc + EI[2])/(EI[3]*Tc2 + EI[4]*Tc + 1.);
        double Fw = (FI[0]*Tc2 + FI[1]*Tc + FI[2])/(FI[3]*Tc2 + FI[4]*Tc + 1.);
        double Eb = Ew + E1*snacl;
        double Fb = Fw + F1*snacl1p5 + F2*snacl + F3*snaclsr;
        double Ib = std::log(std::abs(Eb*(P/P0) + Fb))/Eb;
        double Ib0 = std::log(std::abs(Eb + Fb))/Eb;

        double rhow0 = (DI[0]*Tc2 + DI[1]*Tc + DI[2])/(DI[3]*Tc2 + DI[4]*Tc + 1.);
        double rhob0 = rhow0 + D1*snacl2 + D2*snacl1p5 + D3*snacl + D4*snaclsr;
        double rhobna = rhob0*std::exp(Ib - Ib0);
        double rhobk = rhobna + A2I[0] + A2I[1]*salt + A2I[2]*T + A2I[3]*T2;
        double rhobca = rhobna + B2I[0] + B2I[1]*salt + B2I[2]*T + B2I[3]*T2;

        double Rb = 0.01*(NaCl*rhobna + KCl*rhobk + CaCl2*rhobca);

        return Rb;
    }


    /* Density of brine for an array of (P, T, Sa, NaCl, KCl, CaCl2) values.
     * 
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     *   Sa : double[]
     *     Array of salinity (ppm).
     *   NaCl : double[]
     *     Array of weight percentage of NaCl of total salt in brine.
     *   KCl : double[]
     *     Array of weight percentage of KCl of total salt in brine.
     *   CaCl2 : double[]
     *     Array of weight percentage of CaCl2 of total salt in brine.
     */
    std::vector<double> density(const std::vector<double>& P, const std::vector<double>& T,
                               const std::vector<double>& Sa, const std::vector<double>& NaCl,
                               const std::vector<double>& KCl, const std::vector<double>& CaCl2) {
        if (!are_same_size(P, T, Sa, NaCl, KCl, CaCl2))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = density(P[i], T[i], Sa[i], NaCl[i], KCl[i], CaCl2[i]);
        }

        return out;
    }


    /* Bulk modulus of brine.
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Saturated vapor temperature (C).
     *   Sa : double
     *     Salinity (ppm).
     *   NaCl : double
     *     Weight percentage of NaCl of total salt in brine.
     *   KCl : double
     *     Weight percentage of KCl of total salt in brine.
     *   CaCl2 : double
     *     Weight percentage of CaCl2 of total salt in brine.
     *
     * Returns:
     *   double
     *     Brine bulk density at corresponding (P, T, Sa, NaCl, KCl, CaCl2).
     */
    double bulk_modulus(double P, double T, double Sa, double NaCl, double KCl, double CaCl2) {
        double rhob = density(P, T, Sa, NaCl, KCl, CaCl2);
        double Vb = velocity(P, T, Sa, NaCl, KCl, CaCl2);
        double Kb = rhob*Vb*Vb;

        return Kb;
    }


    /* Bulk modulus of brine for an array of (P, T, Sa, NaCl, KCl, CaCl2) values.
     * 
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     *   Sa : double[]
     *     Array of salinity (ppm).
     *   NaCl : double[]
     *     Array of weight percentage of NaCl of total salt in brine.
     *   KCl : double[]
     *     Array of weight percentage of KCl of total salt in brine.
     *   CaCl2 : double[]
     *     Array of weight percentage of CaCl2 of total salt in brine.
     */
    std::vector<double> bulk_modulus(const std::vector<double>& P, const std::vector<double>& T,
                                    const std::vector<double>& Sa, const std::vector<double>& NaCl,
                                    const std::vector<double>& KCl, const std::vector<double>& CaCl2) {
        if (!are_same_size(P, T, Sa, NaCl, KCl, CaCl2))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = bulk_modulus(P[i], T[i], Sa[i], NaCl[i], KCl[i], CaCl2[i]);
        }

        return out;
    }


    /* Viscosity of brine.
     *
     * Parameters:
     *   T : double
     *     Temperature (C).
     *   Sa : double
     *     Salinity (ppm).
     *
     * Returns:
     *   double
     *     Brine viscosity (centipose) at corresponding (T, Sa).
     */
    double viscosity(double T, double Sa) {
        double salt = Sa*1e-6;
        double eta = 0.1 + 0.333*salt + (1.65 + 91.9*salt*salt*salt)*std::exp(-(0.42*std::pow(std::pow(salt, 0.8)-0.17, 2) + 0.045)*std::pow(T, 0.8));

        return eta;
    }


    /* Viscosity of brine for an array of (T, Sa) values.
     * 
     * Parameters:
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     *   Sa : double[]
     *     Array of salinity (ppm).
     */
    std::vector<double> viscosity(const std::vector<double>& T, const std::vector<double>& Sa) {
        if (!are_same_size(T, Sa))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(T.size());
        for (size_t i=0; i<T.size(); ++i) {
            out[i] = viscosity(T[i], Sa[i]);
        }

        return out;
    }


    /* Resistivity of brine.
     *
     * Limitations:
     *   Sa : Sa > 0
     *
     * Parameters:
     *   T : double
     *     Temperature (C).
     *   Sa : double
     *     Salinity (ppm).
     *
     * Returns:
     *   double
     *     Resistivity of brine at corresponding (T, Sa).
     */
    double resistivity(double T, double Sa) {
        double Rb75 = 0.0123 + 3647.5/powf(Sa, 0.955);
        double Rb = Rb75*(75 + 6.77)/(T + 6.77);

        return Rb;
    }


    /* Resistivity of brine for an array of (T, Sa) values.
     * 
     * Parameters:
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     *   Sa : double[]
     *     Array of salinity (ppm).
     */
    std::vector<double> resistivity(const std::vector<double>& T, const std::vector<double>& Sa) {
        if (!are_same_size(T, Sa))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(T.size());
        for (size_t i=0; i<T.size(); ++i) {
            out[i] = resistivity(T[i], Sa[i]);
        }

        return out;
    }


    /* Solubility of petroleum gas in brine (GWR).
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Temperature (C).
     *   Sa : double
     *     Salinity (ppm).
     *
     *  Returns:
     *    double
     *      Solubility of petroleum gas in brine at corresponding (P, T, Sa).
     */
    double gas_solubility(double P, double T, double Sa) {
        double T2 = T*T;
        double T3 = T2*T;
        double T4 = T3*T;
        double A = 8.15839 - 6.12265e-2*T + 1.91663e-4*T2 - 2.1654e-7*T3;
        double B = 1.01021e-2 - 7.44241e-5*T + 3.05553e-7*T2 - 2.94883e-10*T3;
        double C = (-9.02505 + 0.130237*T - 8.53425e-4*T2 + 2.34122e-6*T3 - 2.37049e-9*T4)*1e-7;

        double Rsw = A + B*P + C*P*P;
        double Rsb = Rsw*std::pow(10, -(0.0840655e-4*Sa*std::pow(T, -0.285854)));

        return Rsb;
    }


    /* Gas-brine for an array of (P, T, Sa) values.
     * 
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     *   Sa : double[]
     *     Array of salinity (ppm).
     */
    std::vector<double> gas_solubility(const std::vector<double>& P, const std::vector<double>& T,
                                       const std::vector<double>& Sa) {
        if (!are_same_size(P, T, Sa))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = gas_solubility(P[i], T[i], Sa[i]);
        }

        return out;
    }


} // END NAMESPACE BRINE
