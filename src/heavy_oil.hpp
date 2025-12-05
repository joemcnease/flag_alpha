#pragma once


#include <vector>


namespace HeavyOil {


    namespace V2011 {

        double velocity_p(double P, double T, double G, double rho0, double Rs);
        double velocity_s(double P, double T, double G, double rho0, double Rs);
        double density(double T, double rho0);
        double viscosity(double T, double rho0);

    };


    namespace VT2011 {

        double velocity_p(double P, double T, double G, double rho0, double Rs);
        std::vector<double> velocity_p(const std::vector<double>& P,
                                       const std::vector<double>& T,
                                       const std::vector<double>& G,
                                       const std::vector<double>& rho0,
                                       const std::vector<double>& Rs);
        double velocity_s(double P, double T, double G, double rho0, double Rs);
        std::vector<double> velocity_s(const std::vector<double>& P,
                                       const std::vector<double>& T,
                                       const std::vector<double>& G,
                                       const std::vector<double>& rho0,
                                       const std::vector<double>& Rs);

    };

    namespace VT2018 {

        double velocity_p(double P, double T, double G, double rho0, double Rs);
        std::vector<double> velocity_p(const std::vector<double>& P,
                                       const std::vector<double>& T,
                                       const std::vector<double>& G,
                                       const std::vector<double>& rho0,
                                       const std::vector<double>& Rs);
        double velocity_s(double P, double T, double G, double rho0, double Rs);
        std::vector<double> velocity_s(const std::vector<double>& P,
                                       const std::vector<double>& T,
                                       const std::vector<double>& G,
                                       const std::vector<double>& rho0,
                                       const std::vector<double>& Rs);

    };

    
    namespace PTD {
        
        double velocity_p(double P, double T, double rho0);
        std::vector<double> velocity_p(const std::vector<double>& P,
                                       const std::vector<double>& T,
                                       const std::vector<double>& rho0);
        double velocity_s(double P, double T, double rho0);
        std::vector<double> velocity_s(const std::vector<double>& P,
                                       const std::vector<double>& T,
                                       const std::vector<double>& rho0);

    };


    namespace BeggsRobinson {
        double viscosity(double T, double rho0);
    };


    namespace DeGhetto {
        double viscosity(double T, double rho0);
        double viscosity_extra_heavy(double T, double rho0);
    };


    namespace FrequencyDependent {
        double velocity_p(double T, double rho0, double f);
        double velocity_s(double T, double rho0, double f);
    };


    double liquid_point(double rho0);
    double glass_point(double rho0);


} // END NAMESPACE HEAVYOIL
