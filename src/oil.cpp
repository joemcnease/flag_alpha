#include <cmath>
#include <stdexcept>
#include <iostream>

#include "license.hpp"
#include "utils.hpp"
#include "oil.hpp"


namespace Oil {


    double velocity(double P, double T, double G, double rho0, double Rs) {
        if ((T < 20) || (T > 170))
            std::runtime_error("Temperature must be between 20 and 170 C");
        if ((Rs < 0) || (Rs > 600))
            std::runtime_error("GOR must be between 0 and 600 L/L");
        // Standing (1962), McCain (1990)
        double Pbp = bubble_point_pressure(T, G, rho0, Rs);
        if ((P < Pbp) || (P > 150))
            std::runtime_error("Pressure must be between bubble point and 150 Mpa.");
        double Mgas = 0.001223*Rs*G;
        double API = to_API(rho0);
        if ((API < 20) || (API > 80))
            std::runtime_error("API must be between 20 and 80.");
        double Ms = 0.5858 + 0.003393*API;
        double Ns = 0.3935 + 0.002701*API;
        double rhoa1 = Ms + Ns*std::log(G);
        double rhop1 = (rho0 + Mgas)/(1. + Mgas/rhoa1);
        double Mg = -1.79117 + 0.631771*G;
        double Ng = 2.638607 - 0.64689*G;
        double rhoa2 = (Mg +  Ng*(rhop1/rho0))*rhoa1;
        double eps = 0.113; // Effective coefficient
        double rhopv = (rho0 + eps*Mgas)/(1. + Mgas/rhoa2);
        double APIvseu = to_API(rhopv);

        double A1 = 3940.7*std::pow(rhopv, 0.32162) - 2289.41;
        double B1 = 3.26313 + 0.00879*APIvseu;
        double C1 = 19.6028*std::log(APIvseu + 307.7138) - 109.04;
        double D1 = 0.994194 - (4.7229e-5 * APIvseu);
        double E1 = 0.05193*std::exp(-2.17305*rhopv);

        double A2 = 3940.7*std::pow(rhopv, 0.32162) - 2289.41;
        double B2 = 4.0525 - 0.0025*APIvseu;
        double C2 = 9.3216*std::log(APIvseu) - 36.3220;
        double D2 = 0.9913 - (1.e-5 * APIvseu);
        double E2 = 0.0009 +  0.025396*rhopv;

        double A, B, C, D, E;

        if ((APIvseu > 0) && (APIvseu < 100)) {
            A = A1;
            B = B1;
            C = C1;
            D = D1;
            E = E1;
        } else if ((APIvseu > 100) && (APIvseu < 200)) {
            A = A1;
            B = 0.01*((200. - APIvseu)*B1 + (APIvseu - 100.)*B2);
            C = 0.01*((200. - APIvseu)*C1 + (APIvseu - 100.)*C2);
            D = 0.01*((200. - APIvseu)*D1 + (APIvseu - 100.)*D2);
            E = 0.01*((200. - APIvseu)*E1 + (APIvseu - 100.)*E2);
        } else if (APIvseu > 200) {
            A = A2;
            B = B2;
            C = C2;
            D = D2;
            E = E2;
        }

        double V = 0.001*(A - B*T + C*(1.-std::pow(D, P))/(1.-D) + E*T*P);

        return V;
    }

    std::vector<double> velocity(const std::vector<double>& P,
                                            const std::vector<double>& T,
                                            const std::vector<double>& G,
                                            const std::vector<double>& rho0,
                                            const std::vector<double>& Rs) {
        if (!are_same_size(P, T, G, rho0, Rs))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = velocity(P[i], T[i], G[i], rho0[i], Rs[i]);
        }

        return out;
    }


    double density(double P, double T, double G, double rho0, double Rs) {
        double API = to_API(rho0);
        if ((API < 20) || (API > 80))
            std::runtime_error("API must be between 20 and 80.");

        // Standing (1962), McCain (1990)
        double Pbp = bubble_point_pressure(T, G, rho0, Rs);
        
        // Katz (1942), McCain (1990)
        double rhoa = (38.52*std::pow(10., -0.00326*API) + (94.75 - 33.93*std::log10(API))*std::log10(G))*0.0160187;

        double rhop0 = rho0;
        double rhoe0 = rho0;
        double Wg = 0;
        if (Rs > 0) {
            double Wg = 0.001223*Rs*G;
            rhop0 = (rho0 + Wg)/(1 + Wg/rhoa);
            double m0 = -1.2214;
            double m1 =  4.6919;
            double m2 = -2.4740;
            double n0 =  0.6093;
            double n1 = -1.1522;
            double n2 =  0.5403;
            double m = m0 + m1*(rhop0/rho0) + m2*std::pow(rhop0/rho0, 2);
            if (m > 1) m = 1;
            double n = n0 + n1*(rhop0/rho0) + n2*std::pow(rhop0/rho0, 2);
            if (n < 0) n = 0;
            rhoe0 = rhop0*(m + n*std::log(P/T));
            if (rhoe0 > rhop0) rhoe0 = rhop0;
        }

        double a = 3.87940e-4 + 3.75885e-2 * std::pow(10., -2.6530*rhoe0);
        double b = 1.00763e-6 + 8.86310e-4 * std::pow(10., -3.7645*rhoe0);
        double drhop = a*P*std::exp(-b/a * P); // ????

        double c0 =  1.6976e-4;
        double c1 =  9.3354e-4;
        double c2 = -1.5383e-6;
        double d0 =  1.2969e1;
        double d1 =  4.0137e-3;
        double d2 = -1.1863e-4;
        double c = c0 + c1*(T - 15.56) + c2*std::pow(T - 15.56, 2);
        double d = d0 + d1*(T - 15.56) + d2*std::pow(T - 15.56, 2);
        double drhot = c*std::exp(-d*drhop);

        std::cout << "G: " << G << ", " << "Wg: " << Wg << ", " << "rhoa: " << rhoa
                  << ", " << "rhop0: " << rhop0 << ", " << "rhoe0: " << rhoe0 << ", " << "drhop: " << drhop << ", "
                  << "drhot: " << drhot << std::endl;

        double rho = rhoe0 + drhop - drhot;

        return rho;
    }

    std::vector<double> density(const std::vector<double>& P,
                                           const std::vector<double>& T,
                                           const std::vector<double>& G,
                                           const std::vector<double>& rho0,
                                           const std::vector<double>& Rs) {
        if (!are_same_size(P, T, G, rho0, Rs))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = density(P[i], T[i], G[i], rho0[i], Rs[i]);
        }

        return out;
    }


    double bulk_modulus(double P, double T, double G, double rho0, double Rs) {
        double V = velocity(P, T, G, rho0, Rs);
        double rho = density(P, T, G, rho0, Rs);
        double K = rho*V*V;

        return K;
    }

    std::vector<double> bulk_modulus(const std::vector<double>& P,
                                                const std::vector<double>& T,
                                                const std::vector<double>& G,
                                                const std::vector<double>& rho0,
                                                const std::vector<double>& Rs) {
        if (!are_same_size(P, T, G, rho0, Rs))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = bulk_modulus(P[i], T[i], G[i], rho0[i], Rs[i]);
        }

        return out;
    }


    double bubble_point_pressure(double T, double G, double rho0, double Rs) {
        double Pbp = std::pow(47.103*Rs/G, 0.83)*std::exp(-(4.072/rho0 - 0.00377*T));

        return Pbp;
    }

    std::vector<double> bubble_point_pressure(const std::vector<double>& T,
                                             const std::vector<double>& G,
                                             const std::vector<double>& rho0,
                                             const std::vector<double>& Rs) {
        if (!are_same_size(T, G, rho0, Rs))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(T.size());
        for (size_t i=0; i<T.size(); ++i) {
            out[i] = bubble_point_pressure(T[i], G[i], rho0[i], Rs[i]);
        }

        return out;
    }


    double viscosity(double P, double T, double G, double rho0, double Rs) {
        double Pbp = bubble_point_pressure(T, G, rho0, Rs);
        double B0 = 0.972 + 3.8e-4 * std::pow(2.4*Rs*std::sqrt(G/rho0) + 20.97 + 17.78, 1.175);
        double rhop = rho0/(B0*(1. +  0.001*Rs));
        double Ly = 5.693 - 2.863/rhop;
        double Y = std::pow(10., Ly);
        double Lnu = 0.505*Y*std::pow(17.8 + T, -1.163);
        double nut = std::pow(10., Lnu) - 1.;
        double Li = 18.6*(0.1*std::log10(nut) + std::pow(std::log10(nut) + 2., -0.1) - 0.985);
        double I = std::pow(10., Li);
        double nu = nut + 0.145*P*I;

        return nu;
    }

    std::vector<double> viscosity(const std::vector<double>& P,
                                 const std::vector<double>& T,
                                 const std::vector<double>& G,
                                 const std::vector<double>& rho0,
                                 const std::vector<double>& Rs) {
        if (!are_same_size(P, T, G, rho0, Rs))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = viscosity(P[i], T[i], G[i], rho0[i], Rs[i]);
        }

        return out;
    }


} // END NAMESPACE OIL
