from utils import close_to, close_to_array

import flag


P     = [  5,  10,  15,  20,  25,  30,  35,  40,  50,  60,  20,  30,  40,  50,    60,    70,    80,    90]
T     = [ 20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,    70,    80,    90,   100]

Vw    = [1.48970, 1.49771, 1.50577, 1.51387, 1.52201, 1.53020, 1.53840, 1.54664, 1.56320,
         1.57987, 1.51387, 1.53019, 1.54664, 1.56320, 1.66332, 1.68457, 1.70283, 1.71844]
rhow  = [1.00042, 1.00269, 1.00493, 1.00715, 1.00934, 1.01151, 1.01365, 1.01577, 1.01995,
         1.02404, 1.00714, 1.01150, 1.01577, 1.01995, 1.00246, 1.00064, 0.99847, 0.99601]
Kw    = [ 2.2201,  2.2492,  2.2785,  2.3082,  2.3381,  2.3684,  2.3990,  2.4298,  2.4923,
          2.5560,  2.3082,  2.3684,  2.4298,  2.4923,  2.7734,  2.8396,  2.8952,  2.9413]

Tnu   = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
nuw   = [0.98080, 0.65316, 0.46391, 0.34609, 0.26966, 0.21870, 0.18403, 0.16007, 0.14330, 0.13144]


def test_h2o_density():
    rho = flag.h2o.density(P[0], T[0])
    assert close_to(rho, rhow[0])

def test_h2o_density_array():
    rho = flag.h2o.density(P, T)
    assert close_to_array(rho, rhow)

def test_h2o_velocity():
    V = flag.h2o.velocity(P[0], T[0])
    assert close_to(V, Vw[0])

def test_h2o_velocity_array():
    V = flag.h2o.velocity(P, T)
    assert close_to_array(V, Vw)

def test_h2o_bulk_modulus():
    K = flag.h2o.bulk_modulus(P[0], T[0])
    assert close_to(K, Kw[0])

def test_h2o_bulk_modulus_array():
    K = flag.h2o.bulk_modulus(P, T)
    assert close_to_array(K, Kw)

def test_h2o_viscosity():
    nu = flag.h2o.viscosity(T[0])
    assert close_to(nu, nuw[0])

def test_h2o_viscosity_array():
    nu = flag.h2o.viscosity(Tnu)
    assert close_to_array(nu, nuw)
