#include <vector>
#include <stdexcept>

#include "license.hpp"
#include "utils.hpp"
#include "h2o.hpp"
#include "h2o_co2.hpp"


namespace H2O_CO2 {


    /* Velocity of H2O+CO2 mixture.
     *
     * Limitations:
     *   Pressure    : 20 MPa - 140 MPa
     *   Temperature : 20 C   - 200 C
     *   Solubility  : Above bubble point
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Temperature (C).
     *   GWR : double
     *     Gas water ratio (L/L).
     *
     * Returns:
     *   double
     *     Velocity of H2O+CO2 mixture at corresponding (P, T, GWR).
     */
    double velocity(double P, double T, double GWR) {
        double Sv = GWR/28.;
        double Tv = 1.4770636e-6*T*T + -6.626486993e-4*T + 2.871802063e-2;
        double Cco2v = Sv*Tv;
        double Vw = H2O::velocity(P, T);
        double V = Vw + Cco2v;

        return V;
    }


    /* Velocity of H2O+CO2 for an array of (P, T, GWR) values.
     * 
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     *   GWR : double[]
     *     Array of gas-water ratio.
     */
    std::vector<double> velocity(const std::vector<double>& P, const std::vector<double>& T,
                                const std::vector<double>& GWR) {
        if (!are_same_size(P, T, GWR))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = velocity(P[i], T[i], GWR[i]);
        }

        return out;
    }


    /* Density of H2O+CO2 mixture.
     *
     * Limitations:
     *   Pressure    : 20 MPa - 140 MPa
     *   Temperature : 20 C   - 200 C
     *   Solubility  : Above bubble point
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Temperature (C).
     *   GWR : double
     *     Gas water ratio (L/L).
     *
     * Returns:
     *   double
     *     Density of H2O+CO2 mixture at corresponding (P, T, GWR).
     */
    double density(double P, double T, double GWR) {
        double Sp = GWR/28.;
        double Tp = -6.20048e-9*T*T + 1.58708e-6*T + -1.363006207e-4  + 1.224855900e-2;
        double Cco2p = Sp*Tp;
        double rhow = H2O::density(P, T);
        double rhowco2 = rhow + Cco2p;

        return rhowco2;
    }


    /* Density of H2O+CO2 for an array of (P, T, GWR) values.
     * 
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     *   GWR : double[]
     *     Array of gas-water ratio.
     */
    std::vector<double> density(const std::vector<double>& P, const std::vector<double>& T,
                               const std::vector<double>& GWR) {
        if (!are_same_size(P, T, GWR))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = density(P[i], T[i], GWR[i]);
        }

        return out;
    }


    /* Bulk modulus of H2O+CO2 mixture.
     *
     * Limitations:
     *   Pressure    : 20 MPa - 140 MPa
     *   Temperature : 20 C   - 200 C
     *   Solubility  : Above bubble point
     *
     * Parameters:
     *   P : double
     *     Pressure (MPa).
     *   T : double
     *     Temperature (C).
     *   GWR : double
     *     Gas water ratio (L/L).
     *
     * Returns:
     *   double
     *     Bulk modulus of H2O+CO2 mixture at corresponding (P, T, GWR).
     */
    double bulk_modulus(double P, double T, double GWR) {
        double rho = density(P, T, GWR); 
        double V = velocity(P, T, GWR); 
        double K = rho*V*V;

        return K;
    }


    /* Bulk modulus of H2O+CO2 for an array of (P, T, GWR) values.
     * 
     * Parameters:
     *   P : double[]
     *     Array of pressure (MPa).
     *   T : double[]
     *     Array of saturated vapor temperature (C).
     *   GWR : double[]
     *     Array of gas-water ratio.
     */
    std::vector<double> bulk_modulus(const std::vector<double>& P, const std::vector<double>& T,
                                    const std::vector<double>& GWR) {
        if (!are_same_size(P, T, GWR))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = bulk_modulus(P[i], T[i], GWR[i]);
        }

        return out;
    }


} // END NAMESPACE H2O_CO2
