#include <stdexcept>

#include "license.hpp"
#include "utils.hpp"
#include "oil.hpp"
#include "oil_co2.hpp"
#include "oil_hc_co2.hpp"


namespace Oil_HC_CO2 {


    double velocity(double P, double T, double Ghc, double Gco2, double rho0, double Rshc, double Rsco2) {
        double Vdeadoil = Oil::velocity(P, T, Ghc, rho0, Rshc);
        double Voilco2 = Oil_CO2::velocity(P, T, Gco2, rho0, Rsco2);
        double dVoilco2 = Vdeadoil - Voilco2;
        double Cgor = Rsco2/600 * (1. + 500./(Rsco2 + 1.));
        double Crho0 = 1.6457*rho0 - 1.3174;
        double Coilhcco2 = Crho0*Cgor;
        double Voilhc = Oil::velocity(P, T, Ghc, rho0, Rshc);

        double V = Voilhc - dVoilco2 + Coilhcco2;

        return V;
    }

    std::vector<double> velocity(const std::vector<double>& P,
                                 const std::vector<double>& T,
                                 const std::vector<double>& Ghc,
                                 const std::vector<double>& Gco2,
                                 const std::vector<double>& rho0,
                                 const std::vector<double>& Rshc,
                                 const std::vector<double>& Rsco2) {
        if (!are_same_size(P, T, Ghc, Gco2, rho0, Rshc, Rsco2))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = velocity(P[i], T[i], Ghc[i], Gco2[i], rho0[i], Rshc[i], Rsco2[i]);
        }

        return out;
    }


    double density(double P, double T, double Ghc, double Gco2, double rho0, double Rshc, double Rsco2) {
        double rho = 0.;
        return rho;
    }

    std::vector<double> density(const std::vector<double>& P,
                                         const std::vector<double>& T,
                                         const std::vector<double>& Ghc,
                                         const std::vector<double>& Gco2,
                                         const std::vector<double>& rho0,
                                         const std::vector<double>& Rshc,
                                         const std::vector<double>& Rsco2) {
        if (!are_same_size(P, T, Ghc, Gco2, rho0, Rshc, Rsco2))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = density(P[i], T[i], Ghc[i], Gco2[i], rho0[i], Rshc[i], Rsco2[i]);
        }

        return out;
    }


    double bulk_modulus(double P, double T, double Ghc, double Gco2, double rho0, double Rshc, double Rsco2) {
        double V = velocity(P, T, Ghc, Gco2, rho0, Rshc, Rsco2);
        double rho = density(P, T, Ghc, Gco2, rho0, Rshc, Rsco2);
        double K = rho*V*V;

        return K;
    }

    std::vector<double> bulk_modulus(const std::vector<double>& P,
                                     const std::vector<double>& T,
                                     const std::vector<double>& Ghc,
                                     const std::vector<double>& Gco2,
                                     const std::vector<double>& rho0,
                                     const std::vector<double>& Rshc,
                                     const std::vector<double>& Rsco2) {
        if (!are_same_size(P, T, Ghc, Gco2, rho0, Rshc, Rsco2))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = bulk_modulus(P[i], T[i], Ghc[i], Gco2[i], rho0[i], Rshc[i], Rsco2[i]);
        }

        return out;
    }


} // END NAMESPACE OIL
