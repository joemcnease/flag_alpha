"""
Physical properties of gas.
"""
from . import _flaglib

from .types import TFVec
from .registry import ModelRegistry

import numpy as np


GAS = ModelRegistry({
    "rpl": { # Hydrocarbon gas model: covers empirical 1999, light, and global model ranges.
        "velocity"    : _flaglib.gas.hydrocarbon_model.velocity,
        "density"     : _flaglib.gas.hydrocarbon_model.density,
        "bulk_modulus": _flaglib.gas.hydrocarbon_model.bulk_modulus,
        "viscosity"   : _flaglib.gas.viscosity
    },
    "global": {
        "velocity"    : _flaglib.gas.global_model.velocity,
        "density"     : _flaglib.gas.global_model.density,
        "bulk_modulus": _flaglib.gas.global_model.bulk_modulus,
    },
    "light": {
        "velocity"    : _flaglib.gas.light_model.velocity,
        "density"     : _flaglib.gas.light_model.density,
        "bulk_modulus": _flaglib.gas.light_model.bulk_modulus,
    },
    "empirical1999": {
        "velocity"    : _flaglib.gas.empirical_model_1999.velocity,
    }
})


def velocity(P: TFVec, T: TFVec, G: TFVec, model: str = "rpl") -> TFVec:
    """Velocity of gas [km/s].

    Application conditions:
        rpl (Hydrocarbon Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   2

        light (Light Gas Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   1

        global (Global Gas Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   2

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    G : TFVec
        Gas gravity
    model : str
        Model used to calculate velocity. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Velocity of gas [km/s].
    """
    v = GAS.call(model, "velocity", P, T, G)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def density(P: TFVec, T: TFVec, G: TFVec, model: str = "rpl") -> TFVec:
    """Density of gas [g/cm^3].

    Application conditions:
        rpl (Hydrocarbon Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   2

        light (Light Gas Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   1

        global (Global Gas Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   2

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    G : TFVec
        Gas gravity
    model : str
        Model used to calculate density. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Density of gas [g/cm^3].
    """
    v = GAS.call(model, "density", P, T, G)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def bulk_modulus(P: TFVec, T: TFVec, G: TFVec, model: str = "rpl") -> TFVec:
    """Bulk modulus of gas [GPa].

    Application conditions:
        rpl (Hydrocarbon Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   2

        light (Light Gas Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   1

        global (Global Gas Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   2

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    G : TFVec
        Gas gravity
    model : str
        Model used to calculate bulk modulus. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Bulk modulus of gas [GPa].
    """
    v = GAS.call(model, "bulk_modulus", P, T, G)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def viscosity(P: TFVec, T: TFVec, G: TFVec, model: str = "rpl") -> TFVec:
    """Viscosity of gas [cP].

    Application conditions:
        rpl (Hydrocarbon Model):
            1. Pressure    [MPa] :  5 <= P <= 110
            2. Temperature [C]   : 20 <= T <= 150
            3. Gas gravity       :  0 <= G <=   2

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    G : TFVec
        Gas gravity
    model : str
        Model used to calculate viscosity. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Viscosity of gas [cP].
    """
    v = GAS.call(model, "viscosity", P, T, G)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v
