from utils import close_to, close_to_array

import flag


P     = [  5,  10,  15,  20,  25,  30,  35,  40,  50,  60,  20,  30,  40,  50,    60,    70,    80,    90]
T     = [ 20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,    70,    80,    90,   100]
Sa    = [2e5, 2e5, 2e5, 2e5, 2e5, 2e5, 2e5, 2e5, 2e5, 2e5, 2e5, 2e5, 2e5, 2e5, 1.5e5, 1.5e5, 1.5e5, 1.5e5]
Na    = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  20,  20,  20,  20,    20,    20,    20,    20]
KCl   = [  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  40,  40,  40,  40,    40,    40,    40,    40]
CaCl2 = [  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  40,  40,  40,  40,    40,    40,    40,    40]

Vb    = [1.72569, 1.73327, 1.74089, 1.74857, 1.75629, 1.76406, 1.77187, 1.77972, 1.79553,
         1.81148, 1.70630, 1.72119, 1.73624, 1.75145, 1.78701, 1.80450, 1.81962, 1.83266]
rhob  = [1.14956, 1.15132, 1.15306, 1.15478, 1.15649, 1.15819, 1.15987, 1.16153, 1.16483,
         1.16807, 1.15258, 1.15598, 1.15933, 1.16263, 1.11040, 1.10961, 1.10865, 1.10755]
Kb    = [ 3.4234,  3.4588,  3.4946,  3.5307,  3.5672,  3.6042,  3.6414,  3.6790,  3.7553,
          3.8330,  3.3557,  3.4246,  3.4948,  3.5664,  3.5459,  3.6131,  3.6708,  3.7199]

Tnu   = [20,  20,  20,  20,    20,    20,    20,    20,    20,  20,  40,  60,  80, 100, 120, 140, 160, 180, 200]
Sanu  = [ 0, 3e4, 6e4, 9e4, 1.2e5, 1.5e5, 1.8e5, 2.1e5, 2.4e5, 1e5, 1e5, 1e5, 1e5, 1e5, 1e5, 1e5, 1e5, 1e5, 1e5]
nub   = [0.98080, 1.06369, 1.11906, 1.17442, 1.24236, 1.33229, 1.45095, 1.60224, 1.78723, 1.19515, 0.86909,
         0.66217, 0.52180, 0.42305, 0.35193, 0.29980, 0.26108, 0.23201, 0.20999]


def test_brine_density():
    rho = flag.brine.density(P[0], T[0], Sa[0], Na[0], KCl[0], CaCl2[0])
    assert close_to(rho, rhob[0])

def test_brine_density_array():
    rho = flag.brine.density(P, T, Sa, Na, KCl, CaCl2)
    assert close_to_array(rho, rhob)

def test_brine_velocity():
    V = flag.brine.velocity(P[0], T[0], Sa[0], Na[0], KCl[0], CaCl2[0])
    assert close_to(V, Vb[0])

def test_brine_velocity_array():
    V = flag.brine.velocity(P, T, Sa, Na, KCl, CaCl2)
    assert close_to_array(V, Vb)

def test_brine_bulk_modulus():
    K = flag.brine.bulk_modulus(P[0], T[0], Sa[0], Na[0], KCl[0], CaCl2[0])
    assert close_to(K, Kb[0])

def test_brine_bulk_modulus_array():
    K = flag.brine.bulk_modulus(P, T, Sa, Na, KCl, CaCl2)
    assert close_to_array(K, Kb)

def test_brine_viscosity():
    nu = flag.brine.viscosity(Tnu[0], Sanu[0])
    assert close_to(nu, nub[0])

def test_brine_viscosity_array():
    nu = flag.brine.viscosity(Tnu, Sanu)
    assert close_to_array(nu, nub)
