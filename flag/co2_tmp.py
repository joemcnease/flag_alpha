from . import _flaglib

from .types import TFVec
from .registry import ModelRegistry


CO2 = ModelRegistry({
    "rpl": {
        "velocity"    : _flaglib.co2.velocity,
        "density"     : _flaglib.co2.density,
        "bulk_modulus": _flaglib.co2.bulk_modulus,
        "pressure"    : _flaglib.co2.pressure
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
    return CO2.call(model, "velocity", P, T)


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
    return CO2.call(model, "bulk_modulus", P, T)


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
        Density of pure CO2 [GPa].
    """
    return CO2.call(model, "pressure", rho, T)
