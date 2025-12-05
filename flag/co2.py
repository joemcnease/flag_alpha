"""
Physical properties of CO2.
"""
from . import _flaglib
from . import co2_eos

from .types import TFVec
from .registry import ModelRegistry

import numpy as np


CO2 = ModelRegistry({
    "rpl": {
        "velocity"    : co2_eos.velocity,
        "density"     : co2_eos.density,
        "bulk_modulus": co2_eos.bulk_modulus,
        "pressure"    : co2_eos.pressure
    }
})


def velocity(P: TFVec, T: TFVec, model: str = "rpl") -> TFVec:
    """Velocity of pure CO2 [km/s].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P <= 800
        2. Temperature [C]   : 0 <= T <= 826

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
        Velocity of pure CO2 [km/s].
    """
    if isinstance(P, (list, np.ndarray)):
        res = []
        for i in range(len(P)):
            Ptmp = P[i]*1e6       # To MPa
            Ttmp = T[i] + 273.15  # To K
            res.append(float(CO2.call(model, "velocity", Ptmp, Ttmp))/1e3)
        if isinstance(P, np.ndarray):
            return np.array(res)
        return res
    return CO2.call(model, "velocity", P, T)/1e3


def density(P: TFVec, T: TFVec, model: str = "rpl") -> TFVec:
    """Density of pure CO2 [g/cm^3].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P <= 800
        2. Temperature [C]   : 0 <= T <= 826

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
        Density of pure CO2 [g/cm^3].
    """
    if isinstance(P, (list, np.ndarray)):
        res = []
        for i in range(len(P)):
            Ptmp = P[i]*1e6       # To MPa
            Ttmp = T[i] + 273.15  # To K
            res.append(float(CO2.call(model, "density", Ptmp, Ttmp))/1e3)
        if isinstance(P, np.ndarray):
            return np.array(res)
        return res
    P *= 1e6
    T += 273.15
    return CO2.call(model, "density", P, T)


def bulk_modulus(P: TFVec, T: TFVec, model: str = "rpl"):
    """Bulk modulus of pure CO2 [GPa].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P <= 800
        2. Temperature [C]   : 0 <= T <= 826

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
        Bulk modulus of pure CO2 [GPa].
    """
    if isinstance(P, (list, np.ndarray)):
        res = []
        for i in range(len(P)):
            Ptmp = P[i]*1e6       # To MPa
            Ttmp = T[i] + 273.15  # To K
            res.append(float(CO2.call(model, "bulk_modulus", Ptmp, Ttmp))/1e9)
        if isinstance(P, np.ndarray):
            return np.array(res)
        return res
    P *= 1e6
    T += 273.15
    return CO2.call(model, "bulk_modulus", P, T)/1e9


def pressure(rho: TFVec, T: TFVec, model: str = "rpl"):
    """Pressure of pure CO2 [MPa].

    Application conditions:
        1. Temperature [C]   : 0 <= T <= 826

    Parameters
    ----------
    rho : TFVec
        Density [g/cm^3]
    T : TFVec
        Temperature [C]
    model : str
        Model used to calculate bulk modulus. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Density of pure CO2 [g/cm^3].
    """
    if isinstance(T, (list, np.ndarray)):
        res = []
        for i in range(len(T)):
            rhotmp = rho[i]*1e3  # To kg/m^3
            Ttmp = T[i] + 273.15  # To K
            res.append(float(CO2.call(model, "pressure", rhotmp, Ttmp))/1e6)
        if isinstance(T, np.ndarray):
            return np.array(res)
        return res
    rho *= 1e3
    T += 273.15
    return CO2.call(model, "pressure", rho, T)/1e6
