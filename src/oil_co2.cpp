#include <cmath>
#include <stdexcept>

#include "license.hpp"
#include "utils.hpp"
#include "oil.hpp"
#include "oil_co2.hpp"


namespace Oil_CO2 {


    double velocity(double P, double T, double G, double rho0, double Rs) {
        // Standing (1962), McCain (1990)
        double Pbp = Oil::bubble_point_pressure(T, G, rho0, Rs);
        double Mco2 = 0.001223*Rs*G;
        double API = to_API(rho0);
        double Ms = 0.564125 + 6.79e-6*API;
        double Ns = 0.132216 + 0.000199874*API;
        double rhoa1 = Ms + Ns*std::log(G);
        double rhop1 = (rho0 + Mco2)/(1. + Mco2/rhoa1);
        double Mg = -1.6125525 + 0.706776*G;
        double Ng = 2.615256 - 0.6753932*G;
        double rhoa2 = (Mg +  Ng*(rhop1/rho0))*rhoa1;
        double eps = 0.113; // Effective coefficient
        double rhopv = (rho0 + eps*Mco2)/(1. + Mco2/rhoa2);
        double APIvseu = to_API(rhopv);

        double A1 = 3940.7*std::pow(rhopv, 0.32162) - 2289.41;
        double B1 = 3.26313 + 0.00879*APIvseu;
        double C1 = 19.6028*std::log(APIvseu + 307.7138) - 109.04;
        double D1 = 0.99221 - (4.71742e-5*APIvseu);
        double E1 = 0.050622*std::exp(-1.60696*rhopv);

        double A2 = 3940.7*std::pow(rhopv, 0.32162) - 2289.41;
        double B2 = 4.0525 - 0.0025*APIvseu;
        double C2 = 9.880896*std::log(APIvseu) - 36.3220;
        double D2 = 0.985352 - (1.482e-5 * APIvseu);
        double E2 = 0.000881 +  0.011597*rhopv;

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

        // Standing (1962), McCain (1990)
        double Pbp = Oil::bubble_point_pressure(T, G, rho0, Rs);

        double Mco2 = 0.001223*Rs*G;

        double a = 3.87940e-4 + 3.75885e-2 * std::pow(10., -2.6530*rho0);
        double b = 1.00763e-6 + 8.86310e-4 * std::pow(10., -3.7645*rho0);
        double drhop = a*P*std::exp(-b/a * P);

        double c0 =  1.6976e4;
        double c1 =  9.3354e-4;
        double c2 = -1.5383e-6;
        double d0 =  1.2969e1;
        double d1 =  4.0137e-3;
        double d2 = -1.1863e-4;
        double c = c0 + c1*(T - 15.56) + c2*std::pow(T - 15.56, 2);
        double d = d0 + d1*(T - 15.56) + d2*std::pow(T - 15.56, 2);
        double drhot = c*std::exp(-d*drhop);

        double rhooil = rho0 + drhop - drhot;

        double rhoeco2 = 0.86476 - 0.001982*T + 0.0074*((1. - std::pow(0.9794, P))/(1. - 0.9794)) + 7.4e-8*T*P;
        double fvco2 = (Mco2/rhoeco2)/(1. + Mco2/rhoeco2);
        double fvoil = 1./(1. + Mco2/rhoeco2);

        double rho = fvco2*rhoeco2 + fvoil*rhooil;

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


} // END NAMESPACE OIL
