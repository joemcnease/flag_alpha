"""
Physical properties of oil.
"""
from . import _flaglib

from .types import TFVec
from .registry import ModelRegistry

import numpy as np


OIL = ModelRegistry({
    "rpl": {
        "velocity"             : _flaglib.oil.velocity,
        "density"              : _flaglib.oil.density,
        "bulk_modulus"         : _flaglib.oil.bulk_modulus,
        "viscosity"            : _flaglib.oil.viscosity,
        "bubble_point_pressure": _flaglib.oil.bubble_point_pressure
    }
})


def velocity(P: TFVec, T: TFVec, G: TFVec, rho0: TFVec, Rs: TFVec, model: str = "rpl") -> TFVec:
    """Velocity of oil [km/s].

    Application conditions:
        1. Pressure           [MPa] : Bubble point <= P    <= 150
        2. Temperature        [C]   :           20 <= T    <= 170
        3. Oil density at STP       :        0.669 <= rho0 <=   0.933
           API gravity              :           20 <= rho0 <=  80
        4. Gas-oil ratio      [L/L] :            0 <= Rs   <= 600

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    G : TFVec
        Gas gravity
    rho0 : TFVec
        Oil density at STP [g/cm^3]
    Rs : TFVec
        Gas-oil ratio [L/L]
    model : str
        Model used to calculate velocity. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Velocity of oil mixture [km/s].
    """
    v = OIL.call(model, "velocity", P, T, G, rho0, Rs)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def density(P, T, G, rho0, Rs, model="rpl"):
    """Density of oil [g/cm^3].

    Application conditions:
        1. Pressure    [MPa] : Bubble point <= P  <= 120
        2. Temperature [C]   :           20 <= T  <= 100
        3. API gravity [L/L] :           20 <= Rs <=  80

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    G : TFVec
        Gas gravity
    rho0 : TFVec
        Oil density at STP [g/cm^3]
    Rs : TFVec
        Gas oil ratio [L/L]
    model : str
        Model used to calculate density. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Density of oil [g/cm^3].
    """
    v = OIL.call(model, "density", P, T, G, rho0, Rs)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def bulk_modulus(P, T, G, rho0, Rs, model="rpl"):
    """Bulk modulus of oil mixture [GPa].

    Application conditions:
        1. Pressure    [MPa] : Bubble point <= P  <= 120
        2. Temperature [C]   :           20 <= T  <= 100
        3. Gas gravity       :            0 <= G  <= 600
        4. API gravity [L/L] :           20 <= Rs <=  80

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    G : TFVec
        Gas gravity
    rho0 : TFVec
        Oil density at STP [g/cm^3]
    Rs : TFVec
        Gas oil ratio [L/L]
    model : str
        Model used to calculate bulk modulus. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Bulk modulus of oil mixture [GPa].
    """
    v = OIL.call(model, "bulk_modulus", P, T, G, rho0, Rs)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def bubble_point_pressure(T, G, rho0, Rs, model="rpl"):
    """Bubble point pressure of oil [MPa].

    Application conditions:

    Parameters
    ----------
    T : TFVec
        Temperature [C]
    G : TFVec
        Gas gravity
    rho0 : TFVec
        Oil density at STP [g/cm^3]
    Rs : TFVec
        Gas oil ratio [L/L]
    model : str
        Model used to calculate bubble point pressure. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Bubble point pressure of oil mixture [MPa].
    """
    v = OIL.call(model, "bubble_point_pressure", T, G, rho0, Rs)
    if isinstance(T, np.ndarray):
        return np.array(v)
    return v


def viscosity(P, T, G, rho0, Rs, model="rpl"):
    """Viscosity of oil [cP].

    Application conditions:
        1. Pressure    [MPa] : Bubble point <= P <=  36
        2. Temperature [C]   :           20 <= T <= 146
        3. Gas gravity       :            5 <= T <= 146

    Parameters
    ----------
    T : TFVec
        Temperature [C]
    Sa : TFVec
        Salinity [ppm]
    model : str
        Model used to calculate viscosity. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Viscosity of oil [cP].
    """
    v = OIL.call(model, "viscosity", P, T, G, rho0, Rs)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v
