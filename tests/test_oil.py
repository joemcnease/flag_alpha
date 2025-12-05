from utils import close_to, close_to_array

import flag


P    = [  40,   40,   40,   40,   40,   40,   40,   60,   60,   60,   60,   60,   60,   60]
T    = [  30,   30,   30,   30,   30,   30,   30,   30,   30,   30,   30,   30,   30,   30]
G    = [0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65]
rho0 = [0.90, 0.88, 0.86, 0.84, 0.82, 0.80, 0.78, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85]
Rs   = [ 100,  100,  100,  100,  100,  100,  100,    0,   50,  100,  150,  200,  250,  300]

Vo   = [1.40406, 1.38973, 1.37507, 1.36005, 1.34467, 1.32890, 1.31272, 1.61914, 1.54397,
        1.46743, 1.39088, 1.30925, 1.24745, 1.21141]
rhoo = [0.818906, 0.800518, 0.781993, 0.763372, 0.744689, 0.725974, 0.707255,
        0.871638, 0.824116, 0.786358, 0.751427, 0.718555, 0.689831, 0.664978]
Ko   = [r*v**2 for (r, v) in zip(rhoo, Vo)]

Pnu    = [ 0.1,   10,   20,   30,   25,   25,   25,   50,   50,   50,   50,   50,   50]
Tnu    = [ 150,  150,  150,  150,   50,   80,  140,  100,  100,  100,  100,  100,  100]
Gnu    = [0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65]
rho0nu = [0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.91, 0.91, 0.91, 0.91, 0.91, 0.91]
Rsnu   = [   0,    0,    0,    0,   50,   50,   50,    0,   30,   60,   90,  120,  150]
nuo    = [1.290224518, 1.490824209, 1.693450160, 1.896076110, 1.528883351,
          0.873656516, 0.468664813, 7.545325146, 2.363260185, 0.929928507,
          0.437899024, 0.277488637, 0.454295531]

rho0p = [0.85, 0.85,  0.85,  0.85,  0.85,  0.85,  0.85,  0.85,  0.85,  0.85,  0.85]
Rsp   = [   0,   50,   100,   150,   200,   250,   100,   100,   100,   100,   100]
Gp    = [0.65, 0.65,  0.65,  0.65,  0.65,  0.65,  0.65,  0.65,  0.65,  0.65,  0.65]
Tp    = [  30,   30,    30,    30,    30,    30,    20,    40,    60,    80,   100]
Pbpo  = [0.00, 8.37, 14.88, 20.83, 26.45, 31.83, 14.33, 15.45, 16.66, 17.96, 19.37]


def test_oil_density():
    rho = flag.oil.density(P[0], T[0], G[0], rho0[0], Rs[0])
    assert close_to(rho, rhoo[0])

def test_oil_density_array():
    rho = flag.oil.density(P, T, G, rho0, Rs)
    assert close_to_array(rho, rhoo)

def test_oil_velocity():
    V = flag.oil.velocity(P[0], T[0], G[0], rho0[0], Rs[0])
    assert close_to(V, Vo[0])

def test_oil_velocity_array():
    V = flag.oil.velocity(P, T, G, rho0, Rs)
    assert close_to_array(V, Vo)

def test_oil_bulk_modulus():
    K = flag.oil.bulk_modulus(P[0], T[0], G[0], rho0[0], Rs[0])
    assert close_to(K, Ko[0])

def test_oil_bulk_modulus_array():
    K = flag.oil.bulk_modulus(P, T, G, rho0, Rs)
    assert close_to_array(K, Ko)

def test_oil_viscosity():
    nu = flag.oil.viscosity(Pnu[0], Tnu[0], Gnu[0], rho0nu[0], Rsnu[0])
    assert close_to(nu, nuo[0])

def test_oil_viscosity_array():
    nu = flag.oil.viscosity(Pnu, Tnu, Gnu, rho0nu, Rsnu)
    assert close_to_array(nu, nuo)

def test_oil_bubble_pressure_point():
    Pbp = flag.oil.bubble_point_pressure(Tp[0], Gp[0], rho0p[0], Rsp[0])
    assert close_to(Pbp, Pbpo[0], tol=1e-1)

def test_oil_bubble_pressure_point_array():
    Pbp = flag.oil.bubble_point_pressure(Tp, Gp, rho0p, Rsp)
    for i in range(len(Pbp)):
        print(abs(Pbp[i]-Pbpo[i]))
    assert close_to_array(Pbp, Pbpo, tol=1e-2)
