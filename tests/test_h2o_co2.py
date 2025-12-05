from utils import close_to, close_to_array

import flag


P   = [ 20,  20,  20,  20,  20,  20,  10,  20,  30,  40,  50,  60,  70,  80, 90, 100]
T   = [ 25,  25,  25,  25,  25,  25,  20,  30,  40,  50,  60,  70,  80,  30, 40,  50]
GWR = [  5,   7,  10,  12,  15,  17,  28,  30,  30,  30,  35,  35,  35,  35, 35,  35]

Vt   = [1.53027, 1.53120, 1.53260, 1.53354, 1.53494, 1.53587, 1.51376, 1.55134,
        1.58268, 1.60959, 1.63157, 1.65029, 1.66602, 1.65224, 1.68274, 1.70876]
rhot = [1.00741, 1.00810, 1.00914, 1.00984, 1.01088, 1.01158, 1.01256, 1.01433,
        1.01459, 1.01419, 1.01476, 1.01333, 1.01149, 1.04039, 1.03963, 1.03841]
Kt   = [r*v*v for r, v in zip(rhot, Vt)]


def test_h2o_co2_density():
    rho = flag.h2o_co2.density(P[0], T[0], GWR[0])
    assert close_to(rho, rhot[0], tol=1e-2)

def test_h2o_co2_density_array():
    rho = flag.h2o_co2.density(P, T, GWR)
    assert close_to_array(rho, rhot, tol=1e-2)

def test_h2o_co2_velocity():
    V = flag.h2o_co2.velocity(P[0], T[0], GWR[0])
    assert close_to(V, Vt[0])

def test_h2o_co2_velocity_array():
    V = flag.h2o_co2.velocity(P, T, GWR)
    assert close_to_array(V, Vt)

def test_h2o_co2_bulk_modulus():
    K = flag.h2o_co2.bulk_modulus(P[0], T[0], GWR[0])
    assert close_to(K, Kt[0], tol=1e-2)

def test_h2o_co2_bulk_modulus_array():
    K = flag.h2o_co2.bulk_modulus(P, T, GWR)
    assert close_to_array(K, Kt, tol=1e-1)
