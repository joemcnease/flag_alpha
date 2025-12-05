#pragma once


#include <vector>


namespace Oil {


    double velocity(double P, double T, double G, double rho0, double Rs);
    std::vector<double> velocity(const std::vector<double>& P,
                                 const std::vector<double>& T,
                                 const std::vector<double>& G,
                                 const std::vector<double>& rho0,
                                 const std::vector<double>& Rs);

    double density(double P, double T, double G, double rho0, double Rs);
    std::vector<double> density(const std::vector<double>& P,
                                const std::vector<double>& T,
                                const std::vector<double>& G,
                                const std::vector<double>& rho0,
                                const std::vector<double>& Rs);

    double bulk_modulus(double P, double T, double G, double rho0, double Rs);
    std::vector<double> bulk_modulus(const std::vector<double>& P,
                                     const std::vector<double>& T,
                                     const std::vector<double>& G,
                                     const std::vector<double>& rho0,
                                     const std::vector<double>& Rs);


    double bubble_point_pressure(double T, double G, double rho0, double Rs);
    std::vector<double> bubble_point_pressure(const std::vector<double>& T,
                                              const std::vector<double>& G,
                                              const std::vector<double>& rho0,
                                              const std::vector<double>& Rs);


    double viscosity(double P, double T, double G, double rho0, double Rs);
    std::vector<double> viscosity(const std::vector<double>& P,
                                  const std::vector<double>& T,
                                  const std::vector<double>& G,
                                  const std::vector<double>& rho0,
                                  const std::vector<double>& Rs);


} // END NAMESPACE OIL
