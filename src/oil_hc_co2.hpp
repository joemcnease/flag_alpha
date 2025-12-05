#pragma once


#include <vector>


namespace Oil_HC_CO2 {


    double velocity(double P, double T, double Ghc, double Gco2, double rho0, double Rshc, double Rsco2);
    std::vector<double> velocity(const std::vector<double>& P,
                                const std::vector<double>& T,
                                const std::vector<double>& Ghc,
                                const std::vector<double>& Gco2,
                                const std::vector<double>& rho0,
                                const std::vector<double>& Rshc,
                                const std::vector<double>& Rsco2);

    double density(double P, double T, double Ghc, double Gco2, double rho0, double Rshc, double Rsco2);
    std::vector<double> density(const std::vector<double>& P,
                                const std::vector<double>& T,
                                const std::vector<double>& Ghc,
                                const std::vector<double>& Gco2,
                                const std::vector<double>& rho0,
                                const std::vector<double>& Rshc,
                                const std::vector<double>& Rsco2);

    double bulk_modulus(double P, double T, double Ghc, double Gco2, double rho0, double Rshc, double Rsco2);
    std::vector<double> bulk_modulus(const std::vector<double>& P,
                                     const std::vector<double>& T,
                                     const std::vector<double>& Ghc,
                                     const std::vector<double>& Gco2,
                                     const std::vector<double>& rho0,
                                     const std::vector<double>& Rshc,
                                     const std::vector<double>& Rsco2);


} // END NAMESPACE OIL_HC_CO2
