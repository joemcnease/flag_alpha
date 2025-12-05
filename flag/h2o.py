"""
Physical properties of H2O.
"""
from . import _flaglib

from .types import TFVec
from .registry import ModelRegistry

import numpy as np


H2O = ModelRegistry({
    "rpl": {
        "velocity"                   : _flaglib.h2o.velocity,
        "density"                    : _flaglib.h2o.density,
        "bulk_modulus"               : _flaglib.h2o.bulk_modulus,
        "viscosity"                  : _flaglib.h2o.viscosity,
        "saturated_vapor_pressure"   : _flaglib.h2o.saturated_vapor_pressure,
        "saturated_vapor_temperature": _flaglib.h2o.saturated_vapor_temperature,
        "gas_solubility"             : _flaglib.h2o.gas_solubility
    }
})


def velocity(P: TFVec, T: TFVec, model: str = "rpl") -> TFVec:
    """Velocity of pure water [km/s].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P <= 200
        2. Temperature [C]   : 0 <= T <= 300

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    model : str
        Model used to calculate velocity. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Velocity of pure water [km/s].
    """
    v = H2O.call(model, "velocity", P, T)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def density(P: TFVec, T: TFVec, model: str = "rpl") -> TFVec:
    """Density of pure water [g/cm^3].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P <= 200
        2. Temperature [C]   : 0 <= T <= 300

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    model : str
        Model used to calculate density. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Density of pure water [g/cm^3].
    """
    v = H2O.call(model, "density", P, T)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def bulk_modulus(P: TFVec, T: TFVec, model: str = "rpl"):
    """Bulk modulus of pure water [GPa].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P <= 200
        2. Temperature [C]   : 0 <= T <= 300

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    model : str
        Model used to calculate bulk modulus. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Bulk modulus of pure water [GPa].
    """
    v = H2O.call(model, "bulk_modulus", P, T)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def viscosity(T: TFVec, model: str = "rpl"):
    """Viscosity of pure water [cP].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P <= 200
        2. Temperature [C]   : 0 <= T <= 300

    Model details:
        From Kestin et al. (1981)

    Parameters
    ----------
    T : TFVec
        Temperature [C]
    model : str
        Model used to calculate viscosity. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Viscosity of pure water [cP].
    """
    v = H2O.call(model, "viscosity", T)
    if isinstance(T, np.ndarray):
        return np.array(v)
    return v


def saturated_vapor_pressure(T: TFVec, model: str = "rpl"):
    """Saturated vapor pressure of pure water [MPa].

    Application conditions:
        1. Temperature [C] : 0.01 <= T <= 373.946

    Model details:
        'rpl':
            Antoine (1888)

    Parameters
    ----------
    T : TFVec
        Temperature [C]
    model : str
        Model used to calculate saturated vapor pressure. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Saturated vapor pressure of pure water [MPa].
    """
    v = H2O.call(model, "saturated_vapor_pressure", T)
    if isinstance(T, np.ndarray):
        return np.array(v)
    return v


def saturated_vapor_temperature(P: TFVec, model: str = "rpl"):
    """Saturated vapor temperature of pure water [C].

    Application conditions:
        1. Pressure [MPa] : 0.0006 <= P <= 21.7175

    Model details:
        'rpl':
            Antoine (1888)

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    model : str
        Model used to calculate saturated vapor temperature. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Saturated vapor temperature of pure water [C].
    """
    v = H2O.call(model, "saturated_vapor_temperature", P)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def gas_solubility(P: TFVec, T: TFVec, model: str = "rpl"):
    """Solubility of petroleum gas in pure water.

    Application conditions:
        1. Pressure    [MPa] : 6.895 <= P <= 68.948
        2. Temperature [C]   :     0 <= T <= Vapor point

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    model : str
        Model used to calculate solubility. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Solubility of petroleum gas in pure water.
    """
    v = H2O.call(model, "gas_solubility", P, T)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v
