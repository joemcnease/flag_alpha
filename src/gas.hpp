#pragma once


#include <vector>


// Global gas model coefficients
static const double ASI_GM[] = { -6.70e-2, 1.67e-2, -1.41e-5,  8.03e-5 };
static const double BSI_GM[] = {  1.67e-2, 8.51e-4,  5.05e-7,  7.45e-8 };
static const double ADI_GM[] = { -3.28e-1, 3.02e-2, -8.27e-7, -5.05e-4 };
static const double BDI_GM[] = {  1.75e-2, 9.11e-4,  2.13e-7,  1.11e-7 };

// Light gas model coefficients
static const double ASI_LM[] = { -2.6243e-1,  2.3800e-2,  2.414e-4, -1.1900e-5 };
static const double BSI_LM[] = {  1.3396e-2,  9.4250e-4,  7.730e-6,  3.7600e-7 };
static const double ADI_LM[] = {  1.50526  , -7.4880e-2, -3.776e-3,  1.9825e-4 };
static const double BDI_LM[] = {  1.7075e-2,  8.1157e-4, -9.950e-7,  6.8700e-7 };

// Hydrocarbon gas density coefficients
static const double ASIJ_HM[5][7] = {
    { -2.0598614e-1,  4.0851307e-1, -2.7072891e-1,  5.1693180e-2,           0., 0., 0. },
    {  4.1998830e-2, -9.2103840e-2,  1.0739269e-1, -4.8468410e-2,  7.540384e-3, 0., 0. },
    {            0.,            0.,            0.,            0.,           0., 0., 0. },
    { -7.5880560e-5,  2.2633320e-4, -2.7244020e-4,  1.2971577e-4, -2.119653e-5, 0., 0. },
    {  8.0300000e-5,            0.,            0.,            0.,           0., 0., 0. },
};
static const double BSIJ_HM[5][7] = {
    {  2.8954460e-2, -4.2402970e-2,  4.6576300e-2, -1.9443810e-2, 2.840336e-3, 0., 0. },
    {  8.5100000e-4,            0.,            0.,            0.,          0., 0., 0. },
    {            0.,            0.,            0.,            0.,          0., 0., 0. },
    {  5.0500000e-7,            0.,            0.,            0.,          0., 0., 0. },
    {  7.4500000e-8,            0.,            0.,            0.,          0., 0., 0. },
};

// Hydrocarbon gas velocity coefficients
static const double ASVIJ_HM1[5][7] = {
	{-2.54972276e+02,  1.93079534e+03, -5.99394395e+03,  9.77725815e+03, -8.85452039e+03,  4.22727918e+03, -8.32154107e+02},
	{ 2.61817807e+00, -2.00437300e+01,  6.34890443e+01, -1.05778613e+02,  9.79384409e+01, -4.78441143e+01,  9.64456106e+00},
	{ 2.95908884e-01, -2.28427742e+00,  7.23550906e+00, -1.20550265e+01,  1.11615238e+01, -5.45253951e+00,  1.09913938e+00},
	{-1.48069766e-02,  1.14319799e-01, -3.62122058e-01,  6.03328800e-01, -5.58610863e-01,  2.72888170e-01, -5.50096215e-02},
	{ 6.63667397e-01, -5.12738305e+00,  1.62411213e+01, -2.70592083e+01,  2.50536154e+01, -1.22389945e+01,  2.46717347e+00},
};
static const double BSVIJ_HM1[5][7] = {
	{-1.35291937e+01,  1.04575376e+02, -3.31207858e+02,  5.51822888e+02, -5.10922501e+02,  2.49591828e+02, -5.03134745e+01},
	{ 4.25753230e-01, -3.28201648e+00,  1.03958739e+01, -1.73204863e+01,  1.60367146e+01, -7.83412926e+00,  1.57922744e+00},
	{ 9.47545847e-03, -7.31460830e-02,  2.31692150e-01, -3.86020525e-01,  3.57409191e-01, -1.74598718e-01,  3.51961368e-02},
	{ 4.30815873e-04, -3.32550630e-03,  1.05336291e-02, -1.75499991e-02,  1.62492162e-02, -7.93793887e-03,  1.60015369e-03},
	{-9.83516946e-03,  7.59848702e-02, -2.40684083e-01,  4.01001917e-01, -3.71280182e-01,  1.81374864e-01, -3.65620929e-02},
};
static const double ADIJ_HM1[5][7] = {
	{ 2.20671437e+03, -1.70370759e+04,  5.39653879e+04, -8.99113214e+04,  8.32472124e+04, -4.06672713e+04,  8.19783136e+03},
	{-1.22540698e+02,  9.46150380e+02, -2.99695632e+03,  4.99320608e+03, -4.62311621e+03,  2.25844825e+03, -4.55264818e+02},
	{-4.62863275e+00,  3.57308680e+01, -1.13178468e+02,  1.88565783e+02, -1.74589535e+02,  8.52891019e+01, -1.71928347e+01},
	{ 2.42962297e-01, -1.87555433e+00,  5.94086788e+00, -9.89803469e+00,  9.16440537e+00, -4.47692299e+00,  9.02471657e-01},
	{ 2.07923544e-01, -1.60742949e+00,  5.09232278e+00, -8.48428017e+00,  7.85543647e+00, -3.83747584e+00,  7.73569968e-01},
};
static const double BDIJ_HM1[5][7] = {
	{-7.97889409e+00,  6.17287615e+01, -1.95494107e+02,  3.25711242e+02, -3.01569953e+02,  1.47320574e+02, -2.96973261e+01},
	{ 2.57267055e-01, -1.98133208e+00,  6.27592170e+00, -1.04562653e+01,  9.68126066e+00, -4.72941306e+00,  9.53369367e-01},
	{-1.48364421e-03,  1.14482123e-02, -3.62578842e-02,  6.04089844e-02, -5.59315499e-02,  2.73232394e-02, -5.50790111e-03},
	{ 1.02225412e-03, -7.89245634e-03,  2.49995640e-02, -4.16515829e-02,  3.85644223e-02, -1.88391873e-02,  3.79766027e-03},
	{-1.58908898e-02,  1.22770350e-01, -3.88878328e-01,  6.47907218e-01, -5.99885186e-01,  2.93051177e-01, -5.90741414e-02},
};
static const double ASVIJ_HM2[5][7] = {
	{ 1.06786966e+00, -3.32124585e+00,  3.00903550e+00, -1.16912827e+00,  1.64695819e-01,  0.00000000e+00,  0.00000000e+00},
	{ 1.11288800e-02,  3.17468200e-02, -2.86932800e-02,  1.11279280e-02, -1.56758700e-03,  0.00000000e+00,  0.00000000e+00},
	{ 2.41400000e-04,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
	{ 6.46872100e-05, -1.89150500e-04,  1.64886400e-04, -6.11909900e-05,  8.08807300e-06,  0.00000000e+00,  0.00000000e+00},
	{ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
};
static const double BSVIJ_HM2[5][7] = {
	{ 6.71099412e-02, -1.46913932e-01,  1.55862460e-01, -6.90108336e-02,  1.03765416e-02,  0.00000000e+00,  0.00000000e+00},
	{-1.13231950e-03,  5.19845300e-03, -4.69826825e-03,  1.82213525e-03, -2.56680450e-04,  0.00000000e+00,  0.00000000e+00},
	{ 7.73000000e-06,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
	{-1.72632880e-06,  5.26738400e-06, -4.76053600e-06,  1.84623520e-06, -2.60082960e-07,  0.00000000e+00,  0.00000000e+00},
	{ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
};
static const double ADIJ_HM2[5][7] = {
	{ 1.50526000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
	{ 5.11677504e-01, -1.45813824e+00,  1.30852800e+00, -5.04264384e-01,  7.01243712e-02,  0.00000000e+00,  0.00000000e+00},
	{-3.77600000e-03,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
	{ 1.98250000e-04,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
	{ 8.03000000e-05,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
};
static const double BDIJ_HM2[5][7] = {
	{ 4.67445200e-02, -8.74086325e-02,  1.01667965e-01, -4.77912175e-02,  7.35864200e-03,  0.00000000e+00,  0.00000000e+00},
	{-5.47006296e-04,  3.52497314e-03, -3.30073635e-03,  1.32083018e-03, -1.91019231e-04,  0.00000000e+00,  0.00000000e+00},
	{-4.97500000e-07,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
	{ 6.87000000e-07,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
	{ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00},
};


namespace Gas {

    
    namespace EmpiricalModel1999 {

        double velocity(double P, double T, double G);
        std::vector<double> velocity(const std::vector<double>& P,
                                     const std::vector<double>& T,
                                     const std::vector<double>& G);

    }


    namespace GlobalModel {

        double velocity(double P, double T, double G);
        std::vector<double> velocity(const std::vector<double>& P,
                                     const std::vector<double>& T,
                                     const std::vector<double>& G);

        double density(double P, double T, double G);
        std::vector<double> density(const std::vector<double>& P,
                                    const std::vector<double>& T,
                                    const std::vector<double>& G);

        double bulk_modulus(double P, double T, double G);
        std::vector<double> bulk_modulus(const std::vector<double>& P,
                                         const std::vector<double>& T,
                                         const std::vector<double>& G);
    }


    namespace LightModel {

        double velocity(double P, double T, double G);
        std::vector<double> velocity(const std::vector<double>& P,
                                     const std::vector<double>& T,
                                     const std::vector<double>& G);

        double density(double P, double T, double G);
        std::vector<double> density(const std::vector<double>& P,
                                    const std::vector<double>& T,
                                    const std::vector<double>& G);

        double bulk_modulus(double P, double T, double G);
        std::vector<double> bulk_modulus(const std::vector<double>& P,
                                         const std::vector<double>& T,
                                         const std::vector<double>& G);

    }


    namespace HydrocarbonModel {

        double velocity(double P, double T, double G);
        std::vector<double> velocity(const std::vector<double>& P,
                                     const std::vector<double>& T,
                                     const std::vector<double>& G);

        double density(double P, double T, double G);
        std::vector<double> density(const std::vector<double>& P,
                                    const std::vector<double>& T,
                                    const std::vector<double>& G);

        double bulk_modulus(double P, double T, double G);
        std::vector<double> bulk_modulus(const std::vector<double>& P,
                                         const std::vector<double>& T,
                                         const std::vector<double>& G);

    }


    namespace CombinedModel {

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

    }


    double viscosity(double P, double T, double G); 
    std::vector<double> viscosity(const std::vector<double>& P,
                                  const std::vector<double>& T,
                                  const std::vector<double>& G);


} // END NAMESPACE GAS
