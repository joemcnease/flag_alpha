#include <cmath>
#include <iostream>
#include <stdexcept>

#include "license.hpp"
#include "utils.hpp"
#include "oil.hpp"
#include "gas.hpp"
#include "oil_gas.hpp"


namespace Oil_Gas {


    double velocity(double P, double T, double G, double rho0, double Rs) {
        double API = to_API(rho0); 
        double Wg = 0.001223*Rs*G;
        // double rhoa = (38.52*std::pow(10., -0.00326*API) + (94.75 - 33.93*std::log(API))*std::log(G))*0.0160187;
        double rhoa = (38.52*std::pow(10., -0.00326*API) + (94.75 - 33.93*std::log10(API))*std::log10(G))*0.0160187;
        double rhop0 = (rho0 + Wg)/(1. + Wg/rhoa);
        double APIseu = to_API(rhop0);
        double eps = 0.113;
        double rhopv = (rho0 + eps*Wg)/(1. + Wg/rhoa);
        double APIvseu = to_API(rhopv);

        double A = 3940.7*std::pow(rhopv, 0.32162) - 2289.41;
        double B = 3.26313 + 0.00879*APIvseu;
        double C = 19.6028*std::log(APIvseu + 307.7138) - 109.04;
        double D = 0.99221 - (5.479e-5 * APIvseu);
        double E = 0.10064*std::exp(-2.6087*rhopv);

        double Voil = 0.001*(A - B*T + C*(1.-std::pow(D, P))/(1.-D) + E*T*P);

        double Molgas = Rs/23.6906;
        double Mwgas = 28.8*G;
        double Wgas = Mwgas*Molgas;
        double Moloil = (1000.*rho0)/(50.*rho0*rho0 + 650.*std::pow(rho0, 7));
        double G0 = (1000.*rho0/Moloil)/28.8;
        double Q = 0.72931 - 0.00436*T - 0.01071*P + 0.00006*T*T + 0.00015*P*P - 0.0001*T*P;
        double Gseu = (1000.*rho0 + Wgas)/((Moloil + Molgas)*28.8);
        eps = std::pow(API/APIseu * Gseu/G0, Q);
        double Gvseu = (eps*1000.*rho0 + Wgas)/((Moloil + Molgas)*28.8);

        double Vgas = Gas::GlobalModel::velocity(P, T, Gvseu);

        double Woil = (110. - APIseu)/(110. - 90.);
        double V = Woil*Voil + (1. - Woil)*Vgas;

        std::cout << Gseu << ", " << Gvseu << std::endl;

        if (APIseu < 90) {
            V = Voil;
        } else if (APIseu > 110) {
            V = Vgas;
        }


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
        double rhooil = Oil::density(P, T, G, rho0, Rs);

        double API = to_API(rho0);
        double Wg = 0.001223*Rs*G;
        double rhoa = (38.52*std::pow(10., -0.00326*API) + (94.75 - 33.93*std::log10(API))*std::log10(G))*0.0160187;
        double rhop0 = (rho0 + Wg)/(1 + Wg/rhoa);
        double APIseu = to_API(rhop0);

        double Molgas = Rs/23.6906;
        double Mwgas = 28.8*G;
        double Wgas = Mwgas*Molgas;
        double Moloil = (1000.*rho0)/(50.*rho0*rho0 + 650.*std::pow(rho0, 7));
        double G0 = (1000.*rho0/Moloil)/28.8;
        double Q = 0.72931 - 0.00436*T - 0.01071*P + 0.00006*T*T + 0.00015*P*P - 0.0001*T*P;
        double Gseu = (1000.*rho0 + Wgas)/(Moloil + Molgas)/28.8;
        double eps = std::pow(API/APIseu * Gseu/G0, Q);

        double rhogas = Gas::GlobalModel::density(P, T, Gseu);

        double Woil = (110. - APIseu)/(110. - 90.);
        double rho = Woil*rhooil + (1. - Woil)*rhogas;
        if (APIseu < 90) {
            rho = rhooil;
        } else if (APIseu > 110) {
            rho = rhogas;
        }

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


} // END NAMESPACE OIL_Gas
