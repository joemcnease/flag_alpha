#include <cmath>
#include <vector>
#include <stdexcept>

#include "license.hpp"
#include "utils.hpp"
#include "co2.hpp"


namespace CO2 {


    /* 
     * Pressure of CO2 from the equation of state. See Span and Wagner (1996).
     */
    double pressure(double rho, double T) {
        double R = 0.1889241;
        double Tc = 304.1282;
        double rhoc = 467.6;
        double delta = rho/rhoc;
        double tau = Tc/T;

        double dm1 = delta - 1;
        double dm12 = dm1*dm1;
        double tm12 = (tau - 1)*(tau - 1);

        double phi1 = 0;
        double phi2 = 0;
        double phi3 = 0;
        double phi4 = 0;

        int k = 0;
        int l = 0;
        int m = 0;
        for (int i=0; i<42; ++i) {
            if (i >= 0 && i<7) {
                phi1 += SWn[i]*SWd[i]*std::pow(delta, SWd[i]-1)*std::pow(tau, SWt[i]);
            }
            else if (i >= 7 && i < 34) {
                phi2 += SWn[i]*expf(-std::pow(delta, SWc[k]))*std::pow(delta, SWd[i]-1)*std::pow(tau, SWt[i])*(SWd[i] - SWc[k]*std::pow(delta, SWc[k]));
                k += 1;
            }
            else if (i >= 34 && i < 39) {
                phi3 += SWn[i]*std::pow(delta, SWd[i])*std::pow(tau, SWt[i])*expf(-SWalpha[l]*(delta - SWepsil[l])*(delta - SWepsil[l]) \
                        - SWbeta[l]*(tau-SWgamma[l])*(tau-SWgamma[l]))*(SWd[i]/delta - 2*SWalpha[l]*(delta-SWepsil[l]));
                l += 1;
            }
            else if (i >= 39 && i < 42) {
                double theta = (1 - tau) + SWA[m]*std::pow(dm12, 1/(2*SWbeta2[m])); 
                double nabla = theta*theta + SWB[m]*std::pow(dm12, SWa[m]);
                double psi = expf(-SWC[m]*dm12 - SWD[m]*tm12);
                double dndd = dm1*(SWA[m]*theta*(2/SWbeta2[m])*std::pow(dm12, 1/(2*SWbeta2[m]) - 1) + 2*SWB[m]*SWa[m]*std::pow(dm12, SWa[m] - 1));
                double dnddbi = SWb[m]*std::pow(nabla, SWb[m]-1)*dndd;
                double dpdd = -2*SWC[m]*dm1*psi;

                phi4 += SWn[i]*(std::pow(nabla, SWb[m])*(psi + delta*dpdd) + dnddbi*delta*psi);
                m += 1;
            }
            else continue;
        }

        double phird = phi1 + phi2 + phi3 + phi4;

        double p = rho*R*T*(1 + delta*phird)*1000;

        return p;
    }


    /*
     * Pressure of CO2 from the equation of state. See Span and Wagner (1996).
     */
    std::vector<double> pressure(const std::vector<double>& rho, const std::vector<double>& T) {
        if (!are_same_size(rho, T))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(rho.size());
        for (size_t i=0; i<rho.size(); ++i) {
            out[i] = pressure(rho[i], T[i]);
        }

        return out;
    }


    /*
     * Function to minimize to find CO2 density from pressure and temperature.
     */
    double density_objective(double x, double P, double T) {
        double e = pressure(x, T) - P;

        return e;
    }


    /*
     * Density of CO2 from the equation of state. See Span and Wagner (1996).
     */
    double density(double P, double T) {
        // 1500 always converges to liquid phase.
        // using a small value (~1e-8) will converge to gaseous phase.
        double rho = co2_secant(1500, P, T);

        return rho;
    }


    /*
     * Density of CO2 from the equation of state. See Span and Wagner (1996).
     */
    std::vector<double> density(const std::vector<double>& P, const std::vector<double>& T) {
        if (!are_same_size(P, T))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = density(P[i], T[i]);
        }

        return out;
    }


    double velocity(double P, double T) {
        double rho = density(P, T);

        double R = 0.1889241;
        double Tc = 304.1282;
        double rhoc = 467.6;
        double delta = rho/rhoc;
        double tau = Tc/T;

        double dm1 = delta - 1;
        double dm12 = dm1*dm1;
        double tm12 = (tau - 1)*(tau - 1);

        double phiott = -SWao[2]/(tau*tau);
        for (int i=3; i<8; ++i) {
            phiott -= SWao[i]*SWthetao[i]*SWthetao[i]*std::exp(-SWthetao[i]*tau)*std::pow(1-expf(-SWthetao[i]*tau), -2);
        }

        double phi11 = 0;
        double phi21 = 0;
        double phi31 = 0;
        double phi41 = 0;

        double phi12 = 0;
        double phi22 = 0;
        double phi32 = 0;
        double phi42 = 0;

        double phi13 = 0;
        double phi23 = 0;
        double phi33 = 0;
        double phi43 = 0;

        double phi14 = 0;
        double phi24 = 0;
        double phi34 = 0;
        double phi44 = 0;

        int k = 0;
        int l = 0;
        int m = 0;
        for (int i=0; i<42; ++i) {
            if (i >= 0 && i<7) {
                phi11 += SWn[i]*SWd[i]*std::pow(delta, SWd[i]-1)*std::pow(tau, SWt[i]);
                phi12 += SWn[i]*SWd[i]*(SWd[i]-1)*std::pow(delta, SWd[i]-2)*std::pow(tau, SWt[i]);
                phi13 += SWn[i]*SWd[i]*SWt[i]*std::pow(delta, SWd[i]-1)*std::pow(tau, SWt[i]-1);
                phi14 += SWn[i]*SWt[i]*(SWt[i]-1)*std::pow(delta, SWd[i])*std::pow(tau, SWt[i]-2);
            }
            else if (i >= 7 && i < 34) {
                phi21 += SWn[i]*expf(-std::pow(delta, SWc[k]))*std::pow(delta, SWd[i]-1)*std::pow(tau, SWt[i])*(SWd[i] - SWc[k]*std::pow(delta, SWc[k]));
                phi22 += SWn[i]*expf(-std::pow(delta, SWc[k]))*((std::pow(delta, SWd[i]-2)*std::pow(tau, SWt[i])*((SWd[i]-SWc[i]*std::pow(delta, SWc[i])))));
                // Resume here to finish the co2 velocity function!
                // phi23 += SWn[i]*std::exp(-std::pow(delta, SWc[i]))*std::pow(delta, SWd[i]-1)*SWt[i]*std::pow(SWt)
                phi24 += 0;
                k += 1;
            }
            else if (i >= 34 && i < 39) {
                phi31 += SWn[i]*std::pow(delta, SWd[i])*std::pow(tau, SWt[i])*expf(-SWalpha[l]*(delta - SWepsil[l])*(delta - SWepsil[l]) \
                        - SWbeta[l]*(tau-SWgamma[l])*(tau-SWgamma[l]))*(SWd[i]/delta - 2*SWalpha[l]*(delta-SWepsil[l]));
                l += 1;
            }
            else if (i >= 39 && i < 42) {
                double theta = (1 - tau) + SWA[m]*std::pow(dm12, 1/(2*SWbeta2[m])); 
                double nabla = theta*theta + SWB[m]*std::pow(dm12, SWa[m]);
                double psi = expf(-SWC[m]*dm12 - SWD[m]*tm12);
                double dndd = dm1*(SWA[m]*theta*(2/SWbeta2[m])*std::pow(dm12, 1/(2*SWbeta2[m]) - 1) + 2*SWB[m]*SWa[m]*std::pow(dm12, SWa[m] - 1));
                double dnddbi = SWb[m]*std::pow(nabla, SWb[m]-1)*dndd;
                double dpdd = -2*SWC[m]*dm1*psi;

                phi41 += SWn[i]*(std::pow(nabla, SWb[m])*(psi + delta*dpdd) + dnddbi*delta*psi);

                double d2pdd2 = (2*SWC[m]*(delta-1)*(delta-1) - 1)*2*SWC[m]*psi;
                double dpdt = -2*SWD[m]*(tau-1)*psi;
                double d2pdt2 = (2*SWD[m]*(tau-1)*(tau-1) - 1)*2*SWD[m]*psi;
                double d2pdddt = 4*SWC[m]*SWD[m]*(delta-1)*(tau-1)*psi;

                double d2ndd2 = (1/dm1)*dndd + dm12*(4*SWB[m]*SWa[m]*(SWa[m]-1)*std::pow(dm12, SWa[m]-2)
                        + 2*SWA[m]*SWA[m]*(1/(SWB[m]*SWB[m]))*std::pow(dm12, 1/(2*SWB[m]) - 1)*std::pow(dm12, 1/(2*SWB[m]) - 1)
                        + SWA[m]*theta*(4/SWbeta2[m])*(1/(2*SWbeta2[m]) - 1)*std::pow(dm12, 1/(2*SWB[m]) - 2));

                double dndtbi = 2*theta*SWb[m]*std::pow(nabla, SWb[m]-1);
                //double d2ndddtbi = -SWA[m]*SWb[m]*(2/SWbeta2[m])*std::pow(nabla, SWb[m]-1)*dm1*std::pow(dm12, 1/(2*SWbeta2) - 1) \
                                  - 2*theta*(SWb[m]-1)*std::pow(nabla, SWb[m]-2)*dndd;
                double d2ndt2bi = 2*SWb[m]*std::pow(nabla, SWb[m]-1) + 4*theta*theta*SWb[m]*(SWb[m]-1)*std::pow(nabla, SWb[m]-2);

                m += 1;
            }
            else continue; // Some error has occurred if we get here...
        }

        double phird = phi11 + phi21 + phi31 + phi41;
        double phirdd = phi12 + phi22 + phi32 + phi42;
        double phirdt = phi13 + phi23 + phi33 + phi43;
        double phirtt = phi14 + phi24 + phi34 + phi44;

        double w = std::sqrt(1000*R*T*(1 + 2*delta*phird + (delta*delta)*phirdd \
                    - std::pow(1 + delta*phird - delta*tau*phirdt, 2)/(tau*tau * (phiott + phirtt))));

        return w;
    }


    /*
     * Velocity of CO2 from the equation of state. See Span and Wagner (1996).
     */
    std::vector<double> velocity(const std::vector<double>& P, const std::vector<double>& T) {
        if (!are_same_size(P, T))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = velocity(P[i], T[i]);
        }

        return out;
    }


    /*
     * Bulk modulus of CO2 from the equation of state. See Span and Wagner (1996).
     */
    double bulk_modulus(double P, double T) {
        double rho = density(P, T);
        double v = velocity(P, T);

        return rho*v*v;
    }


    /*
     * Bulk modulus of CO2 from the equation of state. See Span and Wagner (1996).
     */
    std::vector<double> bulk_modulus(const std::vector<double>& P, const std::vector<double>& T) {
        if (!are_same_size(P, T))
            throw std::runtime_error("Input arrays must have same size.");

        ErrorCode lec = check_license();
        if (lec == LICENSE_ERROR)
            throw std::runtime_error("Licencse error.");

        std::vector<double> out(P.size());
        for (size_t i=0; i<P.size(); ++i) {
            out[i] = bulk_modulus(P[i], T[i]);
        }

        return out;
    }


} // END NAMESPACE CO2
