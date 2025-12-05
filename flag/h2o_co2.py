"""
Physical properties of H2O+CO2 mixture.
"""
from . import _flaglib

from .types import TFVec
from .registry import ModelRegistry

import numpy as np


H2O_CO2 = ModelRegistry({
    "rpl": {
        "velocity"    : _flaglib.h2o_co2.velocity,
        "density"     : _flaglib.h2o_co2.density,
        "bulk_modulus": _flaglib.h2o_co2.bulk_modulus,
    }
})


def velocity(P: TFVec, T: TFVec, Rs: TFVec, model: str = "rpl") -> TFVec:
    """Velocity of H2O+CO2 mixture [km/s].

    Application conditions:
        1. Pressure        [MPa] : 20 <= P  <= 140
        2. Temperature     [C]   : 20 <= T  <= 200
        3. Gas water ratio [L/L] : 20 <= Rs <=  80

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    Rs : TFVec
        Gas water ratio [L/L]
    model : str
        Model used to calculate velocity. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Velocity of h2o+gas mixture [km/s].
    """
    v = H2O_CO2.call(model, "velocity", P, T, Rs)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def density(P: TFVec, T: TFVec, Rs: TFVec, model: str = "rpl") -> TFVec:
    """Density of H2O+CO2 mixture [g/cm^3].

    Application conditions:
        1. Pressure        [MPa] : 20 <= P  <= 140
        2. Temperature     [C]   : 20 <= T  <= 200
        3. Gas water ratio [L/L] : 20 <= Rs <=  80

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    Rs : TFVec
        Gas water ratio [L/L]
    model : str
        Model used to calculate density. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Density of h2o+gas mixture [km/s].
    """
    v = H2O_CO2.call(model, "density", P, T, Rs)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def bulk_modulus(P: TFVec, T: TFVec, Rs: TFVec, model: str = "rpl") -> TFVec:
    """Bulk modulus of H2O+CO2 mixture [GPa].

    Application conditions:
        1. Pressure        [MPa] : 20 <= P  <= 140
        2. Temperature     [C]   : 20 <= T  <= 200
        3. Gas water ratio [L/L] : 20 <= Rs <=  80

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    Rs : TFVec
        Gas water ratio [L/L]
    model : str
        Model used to calculate bulk modulus. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Bulk modulus of h2o+gas mixture [GPa].
    """
    v = H2O_CO2.call(model, "bulk_modulus", P, T, Rs)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v
