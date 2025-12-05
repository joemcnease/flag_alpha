#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <iostream>

#include "license.hpp"
#include "utils.hpp"
#include "gas.hpp"


double trunc9(double x) {
    return std::trunc(x*1e9)/1e9;
}


namespace Gas {

    
    double EmpiricalModel1999::velocity(double P, double T, double G) {
        double G2 = G*G;
        double A =  531.910 + 764.8*std::log(G);
        double B = -4.71620 + 10.941*G - 3.3353*G2;
        double C =  17.9530 + 2.2533*G - 4.3637*G2;
        double D =   1.0010 - 0.0271*G + 0.0126*G2;
        double E =  -0.0392 + 0.0656*G - 0.0190*G2;

        double V = 0.001*(A - B*T + C*(1.-std::pow(D, P))/(1.-D) + E*T*P);

        return V;
    }

    std::vector<double> EmpiricalModel1999::velocity(const std::vector<double>& P,
                                                     const std::vector<double>& T,
                                                     const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = EmpiricalModel1999::velocity(P[i], T[i], G[i]);
        }

        return out;
    }


    double GlobalModel::velocity(double P, double T, double G) {
        double R = 0.008314;
        double Ta = T + 273.15;
        double Mw = 28.8*G;
        double Mw2 = Mw*Mw;

        double as = ASI_GM[0] + ASI_GM[1]*Mw + ASI_GM[2]*Mw*Ta + ASI_GM[3]*Mw2;
        double bs = BSI_GM[0] + BSI_GM[1]*Mw + BSI_GM[2]*Mw*Ta + BSI_GM[3]*Mw2;
        double ad = ADI_GM[0] + ADI_GM[1]*Mw + ADI_GM[2]*Mw*Ta + ADI_GM[3]*Mw2;
        double bd = BDI_GM[0] + BDI_GM[1]*Mw + BDI_GM[2]*Mw*Ta + BDI_GM[3]*Mw2;

        double P1 = -(R*Ta + P*bs)/P;
        double Q = as/P;
        double R1 = -as*bs/P;
        double A = (3.*Q - P1*P1)/3.;
        double B = (2.*P1*P1*P1 - 9.*P1*Q + 27.*R1)/27.;
        double C = std::sqrt(B*B/4. + A*A*A/27.);
        double AA = std::cbrt(-B/2. + C);
        double BB1 = std::cbrt(std::abs(-B/2. - C));
        double BB = BB1*sgn(-B/2. - C);
        double Vm = AA + BB - P1/3.;

        double rho = Mw / (1000.*Vm);
        double K = (Vm*R*Ta/std::pow(Vm-bd, 2.)) - 2.*ad/(Vm*Vm);
        double V = std::sqrt(0.001*K/rho);

        return V;
    }

    std::vector<double> GlobalModel::velocity(const std::vector<double>& P,
                                              const std::vector<double>& T,
                                              const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = GlobalModel::velocity(P[i], T[i], G[i]);
        }

        return out;
    }

    double GlobalModel::density(double P, double T, double G) {
        double R = 0.008314;
        double Ta = T + 273.15;
        double Mw = 28.8*G;
        double Mw2 = Mw*Mw;

        double as = ASI_GM[0] + ASI_GM[1]*Mw + ASI_GM[2]*Mw*Ta + ASI_GM[3]*Mw2;
        double bs = BSI_GM[0] + BSI_GM[1]*Mw + BSI_GM[2]*Mw*Ta + BSI_GM[3]*Mw2;
        double ad = ADI_GM[0] + ADI_GM[1]*Mw + ADI_GM[2]*Mw*Ta + ADI_GM[3]*Mw2;
        double bd = BDI_GM[0] + BDI_GM[1]*Mw + BDI_GM[2]*Mw*Ta + BDI_GM[3]*Mw2;

        double P1 = -(R*Ta + P*bs)/P;
        double Q = as/P;
        double R1 = -as*bs/P;
        double A = (3.*Q - P1*P1)/3.;
        double B = (2.*P1*P1*P1 - 9.*P1*Q + 27.*R1)/27.;
        double C = std::sqrt(B*B/4. + A*A*A/27.);
        double AA = std::cbrt(-B/2. + C);
        double BB1 = std::cbrt(std::abs(-B/2. - C));
        double BB = BB1*sgn(-B/2. - C);
        double Vm = AA + BB - P1/3.;

        double rho = Mw / (1000.*Vm);

        return rho;
    }

    std::vector<double> GlobalModel::density(const std::vector<double>& P,
                                             const std::vector<double>& T,
                                             const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = GlobalModel::density(P[i], T[i], G[i]);
        }

        return out;
    }

    double GlobalModel::bulk_modulus(double P, double T, double G) {
        double R = 0.008314;
        double Ta = T + 273.15;
        double Mw = 28.8*G;
        double Mw2 = Mw*Mw;

        double as = ASI_GM[0] + ASI_GM[1]*Mw + ASI_GM[2]*Mw*Ta + ASI_GM[3]*Mw2;
        double bs = BSI_GM[0] + BSI_GM[1]*Mw + BSI_GM[2]*Mw*Ta + BSI_GM[3]*Mw2;
        double ad = ADI_GM[0] + ADI_GM[1]*Mw + ADI_GM[2]*Mw*Ta + ADI_GM[3]*Mw2;
        double bd = BDI_GM[0] + BDI_GM[1]*Mw + BDI_GM[2]*Mw*Ta + BDI_GM[3]*Mw2;

        double P1 = -(R*Ta + P*bs)/P;
        double Q = as/P;
        double R1 = -as*bs/P;
        double A = (3.*Q - P1*P1)/3.;
        double B = (2.*P1*P1*P1 - 9.*P1*Q + 27.*R1)/27.;
        double C = std::sqrt(B*B/4. + A*A*A/27.);
        double AA = std::cbrt(-B/2. + C);
        double BB1 = std::cbrt(std::abs(-B/2. - C));
        double BB = BB1*sgn(-B/2. - C);
        double Vm = AA + BB - P1/3.;

        double K = (Vm*R*Ta/std::pow(Vm-bd, 2.)) - 2.*ad/(Vm*Vm);

        return K;
    }

    std::vector<double> GlobalModel::bulk_modulus(const std::vector<double>& P,
                                                  const std::vector<double>& T,
                                                  const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = GlobalModel::bulk_modulus(P[i], T[i], G[i]);
        }

        return out;
    }


    double LightModel::velocity(double P, double T, double G) {
        double R = 0.008314;
        double Ta = T + 273.15;
        double Mw = 28.8*G;
        double Mw2 = Mw*Mw;

        double as = ASI_LM[0] + ASI_LM[1]*Mw + (ASI_LM[2] + ASI_LM[3]*Mw)*Ta;
        double bs = BSI_LM[0] + BSI_LM[1]*Mw + (BSI_LM[2] + BSI_LM[3]*Mw)*Ta;
        double ad = ADI_LM[0] + ADI_LM[1]*Mw + (ADI_LM[2] + ADI_LM[3]*Mw)*Ta;
        double bd = BDI_LM[0] + BDI_LM[1]*Mw + (BDI_LM[2] + BDI_LM[3]*Mw)*Ta;

        double P1 = -(R*Ta + P*bs)/P;
        double Q = as/P;
        double R1 = -as*bs/P;
        double A = (3.*Q - P1*P1)/3.;
        double B = (2.*P1*P1*P1 - 9.*P1*Q + 27.*R1)/27.;
        double C = std::sqrt(B*B/4. + A*A*A/27.);
        double AA = std::cbrt(-B/2. + C);
        double BB1 = std::cbrt(std::abs(-B/2. - C));
        double BB = BB1*sgn(-B/2. - C);
        double Vm = AA + BB - P1/3.;

        double rho = Mw / (1000*Vm);
        double K = (Vm*R*Ta/std::pow(Vm-bd, 2.)) - 2.*ad/(Vm*Vm);
        double V = std::pow(0.001*K/rho, 0.5);

        return V;
    }

    std::vector<double> LightModel::velocity(const std::vector<double>& P,
                                             const std::vector<double>& T,
                                             const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = LightModel::velocity(P[i], T[i], G[i]);
        }

        return out;
    }

    double LightModel::density(double P, double T, double G) {
        double R = 0.008314;
        double Ta = T + 273.15;
        double Mw = 28.8*G;
        double Mw2 = Mw*Mw;

        double as = ASI_LM[0] + ASI_LM[1]*Mw + (ASI_LM[2] + ASI_LM[3]*Mw)*Ta;
        double bs = BSI_LM[0] + BSI_LM[1]*Mw + (BSI_LM[2] + BSI_LM[3]*Mw)*Ta;
        double ad = ADI_LM[0] + ADI_LM[1]*Mw + (ADI_LM[2] + ADI_LM[3]*Mw)*Ta;
        double bd = BDI_LM[0] + BDI_LM[1]*Mw + (BDI_LM[2] + BDI_LM[3]*Mw)*Ta;

        double P1 = -(R*Ta + P*bs)/P;
        double Q = as/P;
        double R1 = -as*bs/P;
        double A = (3.*Q - P1*P1)/3.;
        double B = (2.*P1*P1*P1 - 9.*P1*Q + 27.*R1)/27.;
        double C = std::sqrt(B*B/4. + A*A*A/27.);
        double AA = std::cbrt(-B/2. + C);
        double BB1 = std::cbrt(std::abs(-B/2. - C));
        double BB = BB1*sgn(-B/2. - C);
        double Vm = AA + BB - P1/3.;

        double rho = Mw / (1000*Vm);

        return rho;
    }

    std::vector<double> LightModel::density(const std::vector<double>& P,
                                            const std::vector<double>& T,
                                            const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = LightModel::density(P[i], T[i], G[i]);
        }

        return out;
    }

    double LightModel::bulk_modulus(double P, double T, double G) {
        double R = 0.008314;
        double Ta = T + 273.15;
        double Mw = 28.8*G;
        double Mw2 = Mw*Mw;

        double as = ASI_LM[0] + ASI_LM[1]*Mw + (ASI_LM[2] + ASI_LM[3]*Mw)*Ta;
        double bs = BSI_LM[0] + BSI_LM[1]*Mw + (BSI_LM[2] + BSI_LM[3]*Mw)*Ta;
        double ad = ADI_LM[0] + ADI_LM[1]*Mw + (ADI_LM[2] + ADI_LM[3]*Mw)*Ta;
        double bd = BDI_LM[0] + BDI_LM[1]*Mw + (BDI_LM[2] + BDI_LM[3]*Mw)*Ta;

        double P1 = -(R*Ta + P*bs)/P;
        double Q = as/P;
        double R1 = -as*bs/P;
        double A = (3.*Q - P1*P1)/3.;
        double B = (2.*P1*P1*P1 - 9.*P1*Q + 27.*R1)/27.;
        double C = std::sqrt(B*B/4. + A*A*A/27.);
        double AA = std::cbrt(-B/2. + C);
        double BB1 = std::cbrt(std::abs(-B/2. - C));
        double BB = BB1*sgn(-B/2. - C);
        double Vm = AA + BB - P1/3.;

        double K = (Vm*R*Ta/std::pow(Vm-bd, 2.)) - 2.*ad/(Vm*Vm);

        return K;
    }

    std::vector<double> LightModel::bulk_modulus(const std::vector<double>& P,
                                                 const std::vector<double>& T,
                                                 const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = LightModel::bulk_modulus(P[i], T[i], G[i]);
        }

        return out;
    }


    double HydrocarbonModel::velocity(double P, double T, double G) {
        if (G < 0 || G > 2)
            throw std::runtime_error("G must be greater than or equal to 0 and less than or equal to 2.");

        double R = 0.008314;
        double Ta = T + 273.15;
        double Mw = 28.8*G;
        double Mw2 = Mw*Mw;

        double asij[5] = {0};
        double bsij[5] = {0};
        double adij[5] = {0};
        double bdij[5] = {0};

        if (G <= 1) {
            for (int i=0; i<5; ++i) {
                for (int j=0; j<7; ++j) {
                    asij[i] += ASVIJ_HM1[i][j]*std::pow(G, j);
                    bsij[i] += BSVIJ_HM1[i][j]*std::pow(G, j);
                    adij[i] += ADIJ_HM1[i][j]*std::pow(G, j);
                    bdij[i] += BDIJ_HM1[i][j]*std::pow(G, j);
                }
            }
        } else if (G > 1) {
            for (int i=0; i<5; ++i) {
                for (int j=0; j<7; ++j) {
                    asij[i] += ASVIJ_HM2[i][j]*std::pow(G, j);
                    bsij[i] += BSVIJ_HM2[i][j]*std::pow(G, j);
                    adij[i] += ADIJ_HM2[i][j]*std::pow(G, j);
                    bdij[i] += BDIJ_HM2[i][j]*std::pow(G, j);
                }
            }
        }

        double as = asij[0] + asij[1]*Mw + asij[2]*Ta + asij[3]*Mw*Ta + asij[4]*Mw2;
        double bs = bsij[0] + bsij[1]*Mw + bsij[2]*Ta + bsij[3]*Mw*Ta + bsij[4]*Mw2;
        double ad = adij[0] + adij[1]*Mw + adij[2]*Ta + adij[3]*Mw*Ta + adij[4]*Mw2;
        double bd = bdij[0] + bdij[1]*Mw + bdij[2]*Ta + bdij[3]*Mw*Ta + bdij[4]*Mw2;

        double P1 = -(R*Ta + P*bs)/P;
        double Q = as/P;
        double R1 = -as*bs/P;
        double A = (3.*Q - P1*P1)/3.;
        double B = (2.*P1*P1*P1 - 9.*P1*Q + 27.*R1)/27.;
        double C = std::sqrt(B*B/4. + A*A*A/27.);
        double AA = std::cbrt(-B/2. + C);
        double BB1 = std::cbrt(std::abs(-B/2. - C));
        double BB = BB1*sgn(-B/2. - C);
        double Vm = AA + BB - P1/3.;

        double rho = Mw / (1000*Vm);
        double K = (Vm*R*Ta/std::pow(Vm-bd, 2.)) - 2.*ad/(Vm*Vm);
        double V = std::sqrt(K/(1000.*rho));

        return V;
    }

    std::vector<double> HydrocarbonModel::velocity(const std::vector<double>& P,
                                                   const std::vector<double>& T,
                                                   const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = HydrocarbonModel::velocity(P[i], T[i], G[i]);
        }

        return out;
    }

    double HydrocarbonModel::density(double P, double T, double G) {
        if (G < 0 || G > 2)
            throw std::runtime_error("G must be greater than or equal to 0 and less than or equal to 2.");

        double R = 0.008314;
        double Ta = T + 273.15;
        double Mw = 28.8*G;
        double Mw2 = Mw*Mw;

        double asij[5] = {0};
        double bsij[5] = {0};

        for (int i=0; i<5; ++i) {
            for (int j=0; j<7; ++j) {
                asij[i] += ASIJ_HM[i][j]*std::pow(G, j);
                bsij[i] += BSIJ_HM[i][j]*std::pow(G, j);
            }
        }

        double as = asij[0] + asij[1]*Mw + asij[2]*Ta + asij[3]*Mw*Ta + asij[4]*Mw2;
        double bs = bsij[0] + bsij[1]*Mw + bsij[2]*Ta + bsij[3]*Mw*Ta + bsij[4]*Mw2;

        double P1 = -(R*Ta + P*bs)/P;
        double Q = as/P;
        double R1 = -as*bs/P;
        double A = (3.*Q - P1*P1)/3.;
        double B = (2.*P1*P1*P1 - 9.*P1*Q + 27.*R1)/27.;
        double C = std::sqrt(B*B/4. + A*A*A/27.);
        double AA = std::cbrt(-B/2. + C);
        double BB1 = std::cbrt(std::abs(-B/2. - C));
        double BB = BB1*sgn(-B/2. - C);
        double Vm = AA + BB - P1/3.;

        double rho = Mw / (1000*Vm);

        return rho;
    }

    std::vector<double> HydrocarbonModel::density(const std::vector<double>& P,
                                                  const std::vector<double>& T,
                                                  const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = HydrocarbonModel::density(P[i], T[i], G[i]);
        }

        return out;
    }

    double HydrocarbonModel::bulk_modulus(double P, double T, double G) {
        if (G < 0 || G > 2)
            throw std::runtime_error("G must be greater than or equal to 0 and less than or equal to 2.");

        double V = HydrocarbonModel::velocity(P, T, G);
        double rho = HydrocarbonModel::density(P, T, G);
        double K = rho*V*V;

        return K;
    }

    std::vector<double> HydrocarbonModel::bulk_modulus(const std::vector<double>& P,
                                                       const std::vector<double>& T,
                                                       const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = HydrocarbonModel::bulk_modulus(P[i], T[i], G[i]);
        }

        return out;
    }


    double CombinedModel::velocity(double P, double T, double G, double rho0, double Rs) {
        return 0.;
    }

    std::vector<double> CombinedModel::velocity(const std::vector<double>& P,
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
            out[i] = CombinedModel::velocity(P[i], T[i], G[i], rho0[i], Rs[i]);
        }

        return out;
    }

    double CombinedModel::density(double P, double T, double G, double rho0, double Rs) {
        return 0.;
    }

    std::vector<double> CombinedModel::density(const std::vector<double>& P,
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
            out[i] = CombinedModel::density(P[i], T[i], G[i], rho0[i], Rs[i]);
        }

        return out;
    }

    double CombinedModel::bulk_modulus(double P, double T, double G, double rho0, double Rs) {
        return 0.;
    }

    std::vector<double> CombinedModel::bulk_modulus(const std::vector<double>& P,
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
            out[i] = CombinedModel::bulk_modulus(P[i], T[i], G[i], rho0[i], Rs[i]);
        }

        return out;
    }


    double viscosity(double P, double T, double G) {
        double Ppr = P/(4.892 - 0.4048*G);
        double Tabs = T + 273.15;
        double Tpr = Tabs/(94.72 + 170.75*G);
        double mug1 = 0.0001*(Tpr*(28. + 48.*G - 5.*G*G) - 6.47*(1./(G*G)) + 35.*(1/G) + 1.14*G - 15.55);
        double mug2 = (1057. - 8.08*Tpr)/Ppr + (796*std::sqrt(Ppr) - 704.)/(std::pow(Tpr - 1., 0.7)*(Ppr + 1.)) - 3.24*Tpr - 38.;
        double mug = 0.001*Ppr*mug1*mug2;

        return mug;
    }

    std::vector<double> viscosity(const std::vector<double>& P,
                                  const std::vector<double>& T,
                                  const std::vector<double>& G) {
        if (!are_same_size(P, T, G))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = viscosity(P[i], T[i], G[i]);
        }

        return out;
    }


} // END NAMESPACE GAS
