"""
Author: Joe McNease
Date  : Thu Dec 19 2024

Compute density and velocity of CO2 at:
    - temperature (K) [216 <= T <= 1100]
    - pressure (MPa)  [  0 <= P <=  800]

References
----------
.. [SW] Span, R, and Wagner, W. A New Equation of State for Carbon Dioxide Covering
   the Fluid Region from the Triple-Point Temperature to 1100 K at Pressures
   up to 800 MPa. United States: N. p., 1996. Web. doi:10.1063/1.555991.
"""
import numpy as np
from scipy.optimize import root_scalar


# Regression coefficients.
# See Span and Wagner (1996)
# --------------------------------------------------------------------------------
n = np.array([0.38856823203161, 0.29385475942740E1, -0.55867188534934E1,
             -0.76753199592477, 0.31729005580416, 0.54803315897767,
              0.12279411220335, 0.21658961543220E1, 0.15841735109724E1,
             -0.23132705405503, 0.58116916431436E-1, -0.55369137205382,
              0.48946615909422, -0.24275739843501E-1, 0.62494790501678E-1,
             -0.12175860225246, -0.37055685270086, -0.16775879700426E-1,
             -0.11960736637987, -0.45619362508778E-1, 0.35612789270346E-1,
             -0.74427727132052E-2, -0.17395704902432E-2, -0.21810121289527E-1,
              0.24332166559236E-1, -0.37440133423463E-1, 0.14338715756878,
             -0.13491969083286, -0.23151225053480E-1, 0.12363125492901E-1,
              0.21058321972940E-2, -0.33958519026368E-3, 0.55993651771592E-2,
             -0.30335118055646E-3, -0.21365488688320E3, 0.26641569149272E5,
             -0.24027212204557E5, -0.28341603423999E3, 0.21247284400179E3 ,
             -0.66642276540751, 0.72608632349897, 0.55068668612842E-1
             ])
d = np.array([1, 1, 1, 1, 2, 2, 3, 1, 2, 4, 5, 5, 5, 6, 6, 6, 1,
              1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8,
              2, 2, 2, 3, 3])
t = np.array([0.00, 0.75, 1.00, 2.00, 0.75, 2.00, 0.75, 1.50, 1.50, 2.50,
              0.00, 1.50, 2.00, 0.00, 1.00, 2.00, 3.00, 6.00, 3.00, 6.00,
              8.00, 6.00, 0.00, 7.00, 12.0, 16.0, 22.0, 24.0, 16.0, 24.0,
              8.00, 2.00, 28.0, 14.0, 1.00, 0.00, 1.00, 3.00, 3.00])
c = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
              2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6])
alpha = np.array([25, 25, 25, 15, 20])
beta = np.array([325, 300, 300, 275, 275])
gamma = np.array([1.16, 1.19, 1.19, 1.25, 1.22])
epsil = np.array([1.00, 1.00, 1.00, 1.00, 1.00])
a = np.array([3.5, 3.5, 3.0])
b = np.array([0.875, 0.925, 0.875])
beta2 = np.array([0.3, 0.3, 0.3])
A = np.array([0.7, 0.7, 0.7])
B = np.array([0.3, 0.3, 1.0])
C = np.array([10.0, 10.0, 12.5])
D = np.array([275, 275, 275])
ao = np.array([8.37304456, -3.70454304, 2.5, 1.99427042,
               0.62105248, 0.41195293, 1.04028922, 0.08327678])
thetao = np.array([0, 0, 0, 3.15163, 6.11190, 6.77708, 11.32384, 27.08792])
# --------------------------------------------------------------------------------


def pressure(rho, T):
    """
    Compute pressure of C02 for density (rho) and temperature (T)
    given by the equation of state [SW].

    Parameters
    ----------
    rho : float
        Density (kg/m^3)
    T : float
        Temperature (K)

    Returns
    -------
    float
        Pressure (Pa)

    References
    ----------
    .. [SW] Span, R, and Wagner, W. A New Equation of State for Carbon Dioxide Covering
       the Fluid Region from the Triple-Point Temperature to 1100 K at Pressures
       up to 800 MPa. United States: N. p., 1996. Web. doi:10.1063/1.555991.
    """
    # p(\delta, \tau) = \rho*R*T*(1 + \delta*\phi_r^s)
    R = 0.1889241       # (kJ/(kg K))
    Tc = 304.1282       # (K)
    rhoc = 467.6        # (kg/m^3)
    delta = rho/rhoc
    tau = Tc/T

    dm1 = delta - 1
    dm12 = dm1**2
    tm12 = (tau - 1)**2
    theta = (1 - tau) + A*dm12**(1/(2*beta2))
    nabla = theta**2 + B*dm12**a
    psi = np.exp(-C*dm12 - D*tm12)

    dndd = dm1*(A*theta*(2/beta2)*(dm12**(1/(2*beta2) - 1)) + 2*B*a*(dm12**(a - 1)))
    dnddbi = b*(nabla**(b-1))*dndd
    dpdd = -2*C*dm1*psi

    # \phi_{\delta}^r = ...
    phi1 = np.sum(n[:7]*d[:7]*(delta**(d[:7]-1))*(tau**t[:7]))
    phi2 = np.sum(n[7:34]*np.exp(-delta**c)*((delta**(d[7:34]-1))*(tau**t[7:34]))*(d[7:34] - c*delta**c))
    phi3 = np.sum(n[34:39]*(delta**d[34:39])*(tau**t[34:39])*np.exp(-alpha*(delta-epsil)**2 - beta*(tau-gamma)**2)*(d[34:39]/delta - 2*alpha*(delta-epsil)))
    phi4 = np.sum(n[39:42]*((nabla**b)*(psi + delta*dpdd) + dnddbi*delta*psi))
    phird = phi1 + phi2 + phi3 + phi4
    p = rho*R*T*(1 + delta*phird)*1000 # Is off by exactly 3 orders of magnitude?

    return p


def density(P, T):
    """
    Compute density of C02 given temperature and pressure by a root
    finding process using the equation of state [SW].

    Parameters
    ----------
    P : float
        Pressure (Pa)
    T : float
        Temperature (K)

    Returns
    -------
    float
        Density (kg/m^3)

    Notes
    -----
    Since [SW] do not give explicit equations for the density of C02
    at a given temperature and pressure, we must use an implicit
    equation and solve for density by a root finding process. By
    the form of the equation for fixed temperature and varying density,
    we can reliably find a root if we use the Newton-Raphson or secant
    method with an initial guess of 1500 kg/m^3.

    References
    ----------
    .. [SW] Span, R, and Wagner, W. A New Equation of State for Carbon Dioxide Covering
       the Fluid Region from the Triple-Point Temperature to 1100 K at Pressures
       up to 800 MPa. United States: N. p., 1996. Web. doi:10.1063/1.555991.
    """
    # The guess is above liquid phase, so always converges to liquid phase.
    # A lower guess around 0 gives the gas phase.
    sol = root_scalar(lambda x: pressure(x, T) - P, x0=1500, method='secant')
    # sol = root_scalar(lambda x: pressure(x, T) - P, x0=1500, method='newton')

    return sol.root


def velocity(P, T):
    """
    Compute velocity of C02 given temperature and pressure by a root
    finding process using the equation of state [SW].

    Parameters
    ----------
    P : float
        Pressure (Pa)
    T : float
        Temperature (K)

    Returns
    -------
    float
        Velocity (m/s)

    Notes
    -----
    Since [SW] do not give explicit equations for the density of C02
    at a given temperature and pressure, we must use an implicit
    equation and solve for density by a root finding process. By
    the form of the equation for fixed temperature and varying density,
    we can reliably find a root if we use the Newton-Raphson method
    with an initial guess of 2000 kg/m^3.

    References
    ----------
    .. [SW] Span, R, and Wagner, W. A New Equation of State for Carbon Dioxide Covering
       the Fluid Region from the Triple-Point Temperature to 1100 K at Pressures
       up to 800 MPa. United States: N. p., 1996. Web. doi:10.1063/1.555991.
    """
    # First get density from temperature and pressure.
    rho = density(P, T)

    R = 0.1889241       # (kJ/(kg K))
    Tc = 304.1282       # (K)
    rhoc = 467.6        # (kg/m^3)
    delta = rho/rhoc
    tau = Tc/T

    dm1 = delta - 1
    dm12 = dm1**2
    tm12 = (tau - 1)**2
    theta = (1 - tau) + A*dm12**(1/(2*beta2))
    nabla = theta**2 + B*dm12**a
    psi = np.exp(-C*dm12 - D*tm12)

    dndd = dm1*(A*theta*(2/beta2)*(dm12**(1/(2*beta2) - 1)) + 2*B*a*(dm12**(a - 1)))
    dnddbi = b*(nabla**(b-1))*dndd
    dpdd = -2*C*dm1*psi

    d2pdd2 = (2*C*(delta-1)**2 - 1)*2*C*psi
    dpdt = -2*D*(tau-1)*psi
    d2pdt2 = (2*D*(tau-1)**2 - 1)*2*D*psi
    d2pdddt = 4*C*D*(delta-1)*(tau-1)*psi

    d2ndd2 = (1/dm1)*dndd + dm12*(4*B*a*(a-1)*(dm12**(a-2)) + 2*(A**2)*((1/B)**2)*(dm12**(1/(2*B) - 1))**2 + A*theta*(4/beta2)*(1/(2*beta2) - 1)*dm12**(1/(2*B) - 2))
    d2ndd2bi = b*((nabla**(b-1))*d2ndd2 + (b-1)*(nabla**(b-2))*(dndd)**2)

    dndtbi = -2*theta*b*(nabla**(b-1))

    d2ndddtbi = -A*b*(2/beta2)*(nabla**(b-1))*dm1*(dm12**(1/(2*beta2) - 1)) - 2*theta*b*(b-1)*(nabla**(b-2))*dndd
    d2ndt2bi = 2*b*(nabla**(b-1)) + 4*(theta**2)*b*(b-1)*(nabla**(b-2))

    phiott = -ao[2]/tau**2 - np.sum(ao[3:8]*(thetao[3:8]**2)*np.exp(-thetao[3:8]*tau)*(1-np.exp(-thetao[3:8]*tau))**(-2))

    # \phi_{\delta}^r = ...
    phi1 = np.sum(n[:7]*d[:7]*(delta**(d[:7]-1))*(tau**t[:7]))
    phi2 = np.sum(n[7:34]*np.exp(-delta**c)*((delta**(d[7:34]-1))*(tau**t[7:34]))*(d[7:34] - c*delta**c))
    phi3 = np.sum(n[34:39]*(delta**d[34:39])*(tau**t[34:39])*np.exp(-alpha*(delta-epsil)**2 - beta*(tau-gamma)**2)*(d[34:39]/delta - 2*alpha*(delta-epsil)))
    phi4 = np.sum(n[39:42]*((nabla**b)*(psi + delta*dpdd) + dnddbi*delta*psi))
    phird = phi1 + phi2 + phi3 + phi4

    # \phi_{\delta\delta}^r = ...
    phi1 = np.sum(n[:7]*d[:7]*(d[:7]-1)*(delta**(d[:7]-2))*(tau**t[:7]))
    phi2 = np.sum(n[7:34]*np.exp(-delta**c)*((delta**(d[7:34]-2)*(tau**t[7:34])*((d[7:34]-c*delta**c)*(d[7:34]-1-c*delta**c)-(c**2)*(delta**c)))))
    phi3 = np.sum(n[34:39]*(tau**t[34:39])*np.exp(-alpha*(delta-epsil)**2 - beta*(tau-gamma)**2)*(-2*alpha*(delta**d[34:39]) +
                                                                                                  4*(alpha**2)*(delta**d[34:39])*(delta-epsil)**2 -
                            4*d[34:39]*alpha*(delta**(d[34:39]-1))*(delta-epsil) + d[34:39]*(d[34:39]-1)*(delta**(d[34:39]-2))))
    phi4 = np.sum(n[39:42]*((nabla**b)*(2*dpdd + delta*d2pdd2) + 2*dnddbi*(psi + delta*dpdd) + d2ndd2bi*delta*psi))
    phirdd = phi1 + phi2 + phi3 + phi4

    # \phi_{\delta\tau}^r = ...
    phi1 = np.sum(n[:7]*d[:7]*t[:7]*(delta**(d[:7]-1))*(tau**(t[:7]-1)))
    phi2 = np.sum(n[7:34]*np.exp(-delta**c)*(delta**(d[7:34]-1))*(t[7:34])*(tau**(t[7:34]-1))*(d[7:34] - c*delta**c))
    phi3 = np.sum(n[34:39]*(delta**d[34:39])*(tau**t[34:39])*np.exp(-alpha*(delta-epsil)**2 - beta*(tau-gamma)**2)*(d[34:39]/delta - 2*alpha*(delta-epsil))*(t[34:39]/tau - 2*beta*(tau-gamma)))
    phi4 = np.sum(n[39:42]*((nabla**b)*(dpdt + delta*d2pdddt) + delta*dnddbi*dpdt + dndtbi*(psi + delta*dpdd) + d2ndddtbi*delta*psi))
    phirdt = phi1 + phi2 + phi3 + phi4

    # \phi_{\tau\tau}^r
    phi1 = np.sum(n[:7]*t[:7]*(t[:7]-1)*(delta**d[:7])*(tau**(t[:7]-2)))
    phi2 = np.sum(n[7:34]*t[7:34]*(t[7:34]-1)*(delta**d[7:34])*(tau**(t[7:34]-2))*np.exp(-delta**c))
    phi3 = np.sum(n[34:39]*(delta**d[34:39])*(tau**t[34:39])*np.exp(-alpha*(delta-epsil)**2 - beta*(tau-gamma)**2)*((t[34:39]/tau - 2*beta*(tau-gamma))**2 - t[34:39]/(tau**2) - 2*beta))
    phi4 = np.sum(n[39:42]*delta*(d2ndt2bi*psi + 2*dndtbi*dpdt + (nabla**b)*d2pdt2))
    phirtt = phi1 + phi2 + phi3 + phi4

    # Explicit equation for sound speed. Again 3 orders of magnitude off?
    w = np.emath.sqrt(1000*R*T*(1 + 2*delta*phird + (delta**2)*phirdd - (1 + delta*phird - delta*tau*phirdt)**2/(tau**2 * (phiott + phirtt))))

    return w


def bulk_modulus(P, T):
    """
    Compute bulk modulus of C02 given pressure and temperature by a root
    finding process using the equation of state [SW].

    Parameters
    ----------
    P : float
        Pressure (Pa)
    T : float
        Temperature (K)

    Returns
    -------
    float
        Velocity (m/s)

    Notes
    -----
    Since [SW] do not give explicit equations for the density of C02
    at a given temperature and pressure, we must use an implicit
    equation and solve for density by a root finding process. By
    the form of the equation for fixed temperature and varying density,
    we can reliably find a root if we use the Newton-Raphson method
    with an initial guess of 2000 kg/m^3.

    References
    ----------
    .. [SW] Span, R, and Wagner, W. A New Equation of State for Carbon Dioxide Covering
       the Fluid Region from the Triple-Point Temperature to 1100 K at Pressures
       up to 800 MPa. United States: N. p., 1996. Web. doi:10.1063/1.555991.
    """
    # First get density from temperature and pressure.
    rho = density(P, T)
    v = velocity(P, T)

    return rho*v*v
