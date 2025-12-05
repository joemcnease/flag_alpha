from utils import close_to, close_to_array

import flag


Pem = [ 40.,  40.,  40.,  50.,  60.,  70.,  80.,  90., 100.]
Tem = [ 30.,  30.,  30.,  30.,  30.,  20.,  40.,  60.,  80.]
Gem = [0.56, 0.76, 0.96, 1.06, 1.26, 1.46, 1.66, 1.86, 0.72]
Vem = [0.65575, 0.81555, 0.92767, 1.05909, 1.19732, 1.33357, 1.35279, 1.34799, 1.16371]

Pgl   = [ 40.,  40.,  40.,  50.,  60.,  70.,  80.,  90., 100.]
Tgl   = [ 30.,  30.,  30.,  30.,  20.,  30.,  30.,  30.,  30.]
Ggl   = [0.56, 0.76, 0.96, 1.06, 1.26, 1.46, 1.66, 1.86, 0.72]
Vgl   = [0.76088, 0.78908, 0.86173, 0.98853, 1.17039, 1.25078, 1.35811, 1.45495, 1.28328]
rhogl = [0.24202, 0.34135, 0.41938, 0.46720, 0.53745, 0.57705, 0.61598, 0.64828, 0.40587]
Kgl   = [ 140.11,  212.54,  311.42,  456.54,  736.20,  902.77, 1136.15, 1372.34,  668.39]

Pli   = [    40.,     40.,     40.,     50.,     60.,     70.,     80.,     90.,    100.]
Tli   = [    30.,     30.,     30.,     30.,     30.,     20.,     40.,     60.,     80.]
Gli   = [   0.56,    0.76,    0.96,    1.06,    1.26,    1.46,    1.66,    1.86,    0.72]
Vli   = [0.59240, 0.85591, 1.02170, 1.15713, 1.29020, 1.41737, 1.42554, 1.43564, 1.20047]
rholi = [0.23077, 0.34557, 0.42825, 0.47326, 0.53064, 0.58134, 0.60264, 0.61916, 0.37578]
Kli   = [  80.99,  253.16,  447.03,  633.68,  883.31, 1167.88, 1224.30, 1276.12,  541.54]

Phc   = [ 40.,  40.,  40.,  40.,  40.,  40.,  40.,  40.,  20.]
Thc   = [ 30.,  30.,  30.,  30.,  30.,  20.,  40.,  60.,  80.]
Ghc   = [0.56, 0.76, 0.96, 1.06, 1.26, 1.46, 1.66, 1.86, 0.72]
Vhc   = [0.76489, 0.84872, 0.95399, 0.99033, 1.06498, 1.12045, 1.08255, 1.06461, 0.50000]
rhohc = [0.24858, 0.36391, 0.44158, 0.46996, 0.51404, 0.55670, 0.56785, 0.57728, 0.17947]
Khc   = [r*v**2 for (r, v) in zip(rhohc, Vhc)]


def test_gas_empirical1999_velocity():
    V = flag.gas.velocity(Pem[0], Tem[0], Gem[0], model="empirical1999")
    assert close_to(V, Vem[0])

def test_gas_empirical1999_velocity_array():
    V = flag.gas.velocity(Pem, Tem, Gem, model="empirical1999")
    assert close_to_array(V, Vem)

def test_gas_global_velocity():
    V = flag.gas.velocity(Pgl[0], Tgl[0], Ggl[0], model="global")
    assert close_to(V, Vgl[0], tol=1e-2)

def test_gas_global_velocity_array():
    V = flag.gas.velocity(Pgl, Tgl, Ggl, model="global")
    assert close_to_array(V, Vgl, tol=1e-2)

def test_gas_global_density():
    rho = flag.gas.density(Pgl[0], Tgl[0], Ggl[0], model="global")
    assert close_to(rho, rhogl[0], tol=1e-2)

def test_gas_global_density_array():
    rho = flag.gas.density(Pgl, Tgl, Ggl, model="global")
    assert close_to_array(rho, rhogl, tol=1e-2)

def test_gas_global_bulk_modulus():
    K = flag.gas.bulk_modulus(Pgl[0], Tgl[0], Ggl[0], model="global")
    assert close_to(K, Kgl[0], tol=1)

def test_gas_global_bulk_modulus_array():
    K = flag.gas.bulk_modulus(Pgl, Tgl, Ggl, model="global")
    assert close_to_array(K, Kgl, tol=10)

def test_gas_light_velocity():
    V = flag.gas.velocity(Pli[0], Tli[0], Gli[0], model="light")
    assert close_to(V, Vli[0], tol=1e-2)

def test_gas_light_velocity_array():
    V = flag.gas.velocity(Pli, Tli, Gli, model="light")
    assert close_to_array(V, Vli, tol=1e-2)

def test_gas_light_density():
    rho = flag.gas.density(Pli[0], Tli[0], Gli[0], model="light")
    assert close_to(rho, rholi[0], tol=1e-2)

def test_gas_light_density_array():
    rho = flag.gas.density(Pli, Tli, Gli, model="light")
    assert close_to_array(rho, rholi, tol=1e-1)

def test_gas_light_bulk_modulus():
    K = flag.gas.bulk_modulus(Pli[0], Tli[0], Gli[0], model="light")
    assert close_to(K, Kli[0], tol=1e-1)

def test_gas_light_bulk_modulus_array():
    K = flag.gas.bulk_modulus(Pli, Tli, Gli, model="light")
    assert close_to_array(K, Kli, tol=10)

def test_gas_hydrocarbon_velocity():
    V = flag.gas.velocity(Phc[0], Thc[0], Ghc[0], model="rpl")
    assert close_to(V, Vhc[0], tol=1e-2)

def test_gas_hydrocarbon_velocity_array():
    V = flag.gas.velocity(Phc, Thc, Ghc, model="rpl")
    assert close_to_array(V, Vhc, tol=1e-2)

def test_gas_hydrocarbon_density():
    rho = flag.gas.density(Phc[0], Thc[0], Ghc[0], model="rpl")
    assert close_to(rho, rhohc[0], tol=1e-2)

def test_gas_hydrocarbon_density_array():
    rho = flag.gas.density(Phc, Thc, Ghc, model="rpl")
    assert close_to_array(rho, rhohc, tol=1e-2)

def test_gas_hydrocarbon_bulk_modulus():
    K = flag.gas.bulk_modulus(Phc[0], Thc[0], Ghc[0], model="rpl")
    assert close_to(K, Khc[0], tol=1e-2)

def test_gas_hydrocarbon_bulk_modulus_array():
    K = flag.gas.bulk_modulus(Phc, Thc, Ghc, model="rpl")
    assert close_to_array(K, Khc, tol=1e-2)
