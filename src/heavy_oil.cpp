#include <cmath>
#include <stdexcept>
#include <iostream>
#include <limits>

#include "license.hpp"
#include "utils.hpp"
#include "oil.hpp"
#include "heavy_oil.hpp"


namespace HeavyOil {


    double VT2018::velocity_p(double P, double T, double G, double rho0, double Rs) {
        double API = to_API(rho0);
        double A = 0.22598*std::pow(Rs, 0.94335) + 0.03086;
        double B = Rs*(2.*API + 40.)/24.;
        double C = -1.56983*std::pow(B, 0.28143) + 0.00821*B + 5.11428;
        double Rsvho = Rs*(1. + A);
        double rhovseuho = Oil::density(P, T, G, rho0, Rsvho);
        double APIvseuho = to_API(rhovseuho);
        double Rsc = Rs*(1.15 + C);

        double SPT = -0.008071 + 0.013442*rhovseuho - 0.0060654*rhovseuho*rhovseuho;
        double dVpnon = 5669749940.5783*std::pow(P, 1.17528)/(APIvseuho*(150. + std::pow(150. + T, 5.32334)));
        double t0PT = -375.59 + 366.74*rhovseuho;
        double dTp = T - t0PT;
        double CPT = -0.0934 + 0.0316*rhovseuho;
        double APT = -0.0576 + 1.211*rhovseuho - 0.528*rhovseuho*rhovseuho;

        double Vplin = SPT*(dTp - std::abs(dTp));
        double Vpnon = APT*(std::exp(CPT*dTp)/(std::exp(CPT*dTp + 1.))) + dVpnon;
        double Vpcho = Oil::velocity(P, T, G, rho0, Rsc);

        double V = Vpcho + Vpnon + Vplin;

        return V;
    }


    double V2011::density(double T, double rho0) {
        return rho0 - 0.00055*(T - 15.56);
    }


    double V2011::velocity_p(double P, double T, double G, double rho0, double Rs) {
        double API = to_API(rho0);
        double A = 3940.70*std::pow(rho0, 0.32162) - 2289.41;
        double B = 3.26313 + 0.00879*API;
        double E = 0.10064*std::exp(-2.6078*rho0);
        double Vpflag = 0.001*(A - B*T + 0.1*E*T);

        double dVp = Vpflag - 1.6820;
        double Vp = Vpflag*(1. + 0.38184*std::exp(18.044*dVp)/(std::exp(18.044*dVp) + 1.));

        double Vpt = 0.038647*std::pow(Vp, 5) - 0.241712*std::pow(Vp, 4) + 0.411407*std::pow(Vp, 3)
                   + 0.19117*Vp*Vp - 0.067414*Vp + 0.712929;

        return Vpt;
    }


    double V2011::velocity_s(double P, double T, double G, double rho0, double Rs) {
        double API = to_API(rho0);
        double A = 3940.70*std::pow(rho0, 0.32162) - 2289.41;
        double B = 3.26313 + 0.00879*API;
        double E = 0.10064*std::exp(-2.6078*rho0);
        double Vpflag = 0.001*(A - B*T + 0.1*E*T);

        double dVp = Vpflag - 1.6820;
        double Vp = Vpflag*(1. + 0.38184*std::exp(18.044*dVp)/(std::exp(18.044*dVp) + 1.));
        double dVs = Vpflag - 1.6281;
        double R = 0.44034*std::exp(16.4651*dVs)/(std::exp(16.4651*dVs) + 1.);

        double Vs = R*Vp;

        return Vs;
    }


    double V2011::viscosity(double T, double rho0) {
        double Tk = T + 273.15;
        double Tglass = glass_point(rho0);
        double Tliquid = liquid_point(rho0);
        double Tg = Tglass + 273.15;
        double Tl = Tliquid + 273.15;
        double Tw = 200. + 273.15;
        double b1 = std::log10(std::log10(10.e15 + 1.));
        double b2 = std::log10(std::log10(10.e3 + 1.));
        double b3 = std::log10(std::log10(2));
        double x = std::log10(Tg);
        double y = std::log10(Tl);
        double z = std::log10(Tw);
        double x2 = x*x;
        double y2 = y*y;
        double z2 = z*z;
        double detA = (y*z2 - z*y2) - x*(z2 - y2) + x2*(z - y);
        if (std::abs(detA) < 1e-10) {
            throw std::runtime_error("Determinant of A is 0. Cannot find inverse.");
        }
        double detA1 = b1*(y*z2 - z*y2) - x*(b2*z2 - b3*y2) + x2*(b2*z - b3*y);
        double detA2 = (b2*z2 - b3*y2) - b1*(z2 - y2) + x2*(b3 - b2);
        double detA3 = (b3*y - b2*z) - x*(b3 - b2) + b1*(z - y);
        double Ae = detA1/detA;
        double Ce = detA2/detA;
        double De = detA3/detA;
        double Y = Ae + Ce*std::log10(Tk) + De*std::pow(std::log10(Tk), 2);
        double X = std::pow(10, Y); 

        double eta = std::pow(10, X) - 1.;

        return eta;
    }


    double VT2011::velocity_p(double P, double T, double G, double rho0, double Rs) {
        double API = to_API(rho0);
        double A = 3940.70*std::pow(rho0, 0.32162) - 2289.41;
        double B = 3.26313 + 0.00879*API;
        double E = 0.10064*std::exp(-2.6078*rho0);
        double Vpflag = 0.001*(A - B*T + 0.1*E*T);

        double SPT = -0.008071 + 0.013442*rho0 - 0.0060654*rho0*rho0;
        double t0PT = -375.59 + 366.74*rho0;
        double dTp = T - t0PT;
        double CPT = -0.0934 + 0.0316*rho0;
        double APT = -0.0576 + 1.211*rho0 - 0.528*rho0*rho0;

        double Vplin = SPT*(dTp - std::abs(dTp));
        double Vpnon = APT*std::exp(CPT*dTp)/(std::exp(CPT*dTp) + 1.);
        double Vp = Vpflag + Vpnon + Vplin;

        double Vpt = 0.038647*std::pow(Vp, 5) - 0.241712*std::pow(Vp, 4) + 0.411407*std::pow(Vp, 3)
                   + 0.19117*Vp*Vp - 0.067414*Vp + 0.712929;

        return Vp;
    }


    double VT2011::velocity_s(double P, double T, double G, double rho0, double Rs) {
        double SST = -0.0116 + 0.0197*rho0 - 0.0092*rho0*rho0;
        double toST = -372.57 + 371.72*rho0;
        double dTs = T - toST;
        double CST = -0.0798 + 0.0254*rho0;
        double AST = -0.2870 + 2.4132*rho0 - 1.1324*rho0*rho0;
        double Vsnon = AST*std::exp(CST*dTs)/(std::exp(CST*dTs) + 1.);
        double Vslin = SST*(dTs - std::abs(dTs));

        double Vs = Vsnon + Vslin;

        return Vs;
    }


    double  PTD::velocity_p(double P, double T, double rho0) {
        double API = to_API(rho0);
        double SPT = -0.008071 + 0.013442*rho0 - 0.0060654*rho0*rho0;
        double dVpnon = 5669749940.5783*std::pow(P, 1.17528)/(API*(150. + std::pow(150. + T, 5.32334)));
        double t0PT = -375.59 + 366.74*rho0;
        double dTp = T - t0PT;
        double CPT = -0.0934 + 0.0316*rho0;
        double APT = 0.0576 + 1.211*rho0 - 0.528*rho0*rho0;

        double Vplin = SPT*(dTp - std::abs(dTp));
        double Vpnon = APT*(std::exp(CPT*dTp)/(std::exp(CPT*dTp) + 1.));
        double VpCTP = Oil::velocity(P, T, 1.0, rho0, 0.0);
        
        double Vp = VpCTP + Vpnon + Vplin;

        return Vp;
    }


    double  PTD::velocity_s(double P, double T, double rho0) {
        double SST = -0.0116 + 0.0197*rho0 - 0.0092*rho0*rho0;
        double t0ST = -372.57 + 371.72*rho0;
        double dTs = T - t0ST;
        double t0PT = -375.59 + 366.74*rho0;
        double dTp = T - t0PT;
        double APT = 0.0576 + 1.211*rho0 - 0.528*rho0*rho0;
        double CPT = -0.0934 + 0.0316*rho0;

        double Vslin = SST*(dTs - std::abs(dTs));
        double Vpnon = APT*(std::exp(CPT*dTp)/(std::exp(CPT*dTp) + 1.));
        double Vsnon = -0.539948 + std::sqrt(0.1168288 + 1.2416*Vpnon)/0.6208;

        double Vs = Vsnon + Vslin;

        return Vs;
    }


    // FrequencyDependent::velocity_p(double T, double rho0) {
    //     double Vpfre = Vpflag + Vpnonf + Vplin;
    // }


    double BeggsRobinson::viscosity(double T, double rho0) {
        double Y = std::pow(10., 5.6926 - 2.863/rho0);
        double X = 0.505*Y*std::pow(17.8 + T, -1.163);
        
        double eta = std::pow(10, X) - 1.;
        if ((1.8*T + 32.) <= 0) {
            eta = std::numeric_limits<double>::infinity();
        }
        return eta;
    }


    double DeGhetto::viscosity(double T, double rho0) {
        double API = to_API(rho0);
        double Tlog = std::log10(1.8*T + 32.);
        double X = std::pow(10, 2.06492 - 0.0179*API - 0.70226*Tlog);

        double eta = std::pow(10, X) - 1.;
        if ((1.8*T + 32.) <= 0) {
            eta = std::numeric_limits<double>::infinity();
        }
        return eta;
    }


    double DeGhetto::viscosity_extra_heavy(double T, double rho0) {
        double API = to_API(rho0);
        double Tlog = std::log10(1.8*T + 32.);
        double X = std::pow(10, 1.90296 - 0.012619*API - 0.61748*Tlog);

        double eta = std::pow(10, X) - 1.;
        if ((1.8*T + 32.) <= 0) {
            eta = std::numeric_limits<double>::infinity();
        }
        return eta;
    }


    double liquid_point(double rho0) {
        return -354.17 + 393.14*rho0;
    }


    double glass_point(double rho0) {
        return -378.18 + 335.31*rho0;
    }


} // END NAMESPACE HEAVYOIL
