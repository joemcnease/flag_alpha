#pragma once


#include <vector>


namespace H2O_CO2 {


    double velocity(double P, double T, double GWR);
    std::vector<double> velocity(const std::vector<double>& P, const std::vector<double>& T,
                                const std::vector<double>& GWR);

    double density(double P, double T, double GWR);
    std::vector<double> density(const std::vector<double>& P, const std::vector<double>& T,
                               const std::vector<double>& GWR);

    double bulk_modulus(double P, double T, double GWR);
    std::vector<double> bulk_modulus(const std::vector<double>& P, const std::vector<double>& T,
                                    const std::vector<double>& GWR);


} // END NAMESPACE H2O_CO2
