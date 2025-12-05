"""
Physical properties of brine.
"""
from . import _flaglib

from .types import TFVec
from .registry import ModelRegistry

import numpy as np


BRINE = ModelRegistry({
    "rpl": {
        "velocity"      : _flaglib.brine.velocity,
        "density"       : _flaglib.brine.density,
        "bulk_modulus"  : _flaglib.brine.bulk_modulus,
        "viscosity"     : _flaglib.brine.viscosity,
        "resistivity"   : _flaglib.brine.resistivity,
        "gas_solubility": _flaglib.brine.gas_solubility
    }
})


def velocity(P: TFVec, T: TFVec, Sa: TFVec, NaCl: TFVec, KCl: TFVec,
             CaCl2: TFVec, model: str = "rpl") -> TFVec:
    """Velocity of brine [km/s].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P  <= 140
        2. Temperature [C]   : 0 <= T  <= 300
        3. Salinity    [ppm] : 0 <= Sa <= 300,000

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    Sa : TFVec
        Salinity [ppm]
    NaCl : TFVec
        Weight percent of NaCl of total salt of brine [%].
    KCl : TFVec
        Weight percent of KCl of total salt of brine [%].
    CaCl2 : TFVec
        Weight percent of CaCl2 of total salt of brine [%].
    model : str
        Model used to calculate velocity. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Velocity of brine [km/s].
    """
    v = BRINE.call(model, "velocity", P, T, Sa, NaCl, KCl, CaCl2)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def density(P: TFVec, T: TFVec, Sa: TFVec, NaCl: TFVec, KCl: TFVec,
            CaCl2: TFVec, model: str = "rpl") -> TFVec:
    """Density of brine [g/cm^3].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P  <= 140
        2. Temperature [C]   : 0 <= T  <= 300
        3. Salinity    [ppm] : 0 <= Sa <= 300,000

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    Sa : TFVec
        Salinity [ppm]
    NaCl : TFVec
        Weight percent of NaCl of total salt of brine [%].
    KCl : TFVec
        Weight percent of KCl of total salt of brine [%].
    CaCl2 : TFVec
        Weight percent of CaCl2 of total salt of brine [%].
    model : str
        Model used to calculate density. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Density of pure water [km/s].
    """
    v = BRINE.call(model, "density", P, T, Sa, NaCl, KCl, CaCl2)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def bulk_modulus(P: TFVec, T: TFVec, Sa: TFVec, NaCl: TFVec, KCl: TFVec,
                 CaCl2: TFVec, model: str = "rpl") -> TFVec:
    """Bulk modulus of brine [GPa].

    Application conditions:
        1. Pressure    [MPa] : 0 <= P  <= 140
        2. Temperature [C]   : 0 <= T  <= 300
        3. Salinity    [ppm] : 0 <= Sa <= 300,000

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    Sa : TFVec
        Salinity [ppm]
    NaCl : TFVec
        Weight percent of NaCl of total salt of brine [%].
    KCl : TFVec
        Weight percent of KCl of total salt of brine [%].
    CaCl2 : TFVec
        Weight percent of CaCl2 of total salt of brine [%].
    model : str
        Model used to calculate bulk modulus. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Bulk modulus of pure water [km/s].
    """
    v = BRINE.call(model, "bulk_modulus", P, T, Sa, NaCl, KCl, CaCl2)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v


def viscosity(T: TFVec, Sa: TFVec, model: str = "rpl") -> TFVec:
    """Viscosity of brine [cP].

    Application conditions:
        1. Temperature [C]   : 0 <= T  <= 300
        2. Salinity    [ppm] : 0 <= Sa <= 300,000

    Model details:
        From Kestin et al. (1981)

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
        Viscosity of pure water [cP].
    """
    v = BRINE.call(model, "viscosity", T, Sa)
    if isinstance(T, np.ndarray):
        return np.array(v)
    return v


def resistivity(T: TFVec, Sa: TFVec, model: str = "rpl") -> TFVec:
    """Resistivity of brine [cP].

    Application conditions:
        1. Temperature [C]   : 0 <= T  <= 300
        2. Salinity    [ppm] : 0 <= Sa <= 300,000

    Model details:
        From Western Atlas. (1995)

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
        Viscosity of pure water [cP].
    """
    v = BRINE.call(model, "viscosity", T, Sa)
    if isinstance(T, np.ndarray):
        return np.array(v)
    return v


def gas_solubility(P: TFVec, T: TFVec, Sa: TFVec, model: str = "rpl") -> TFVec:
    """Solubility of petroleum gas in brine.

    Application conditions:
        1. Pressure    [MPa] : 6.895 <= P  <= 68.948
        2. Temperature [C]   :     0 <= T  <= Vapor point
        3. Salinity    [ppm] :     0 <= Sa <= 300,000

    Parameters
    ----------
    P : TFVec
        Pressure [MPa]
    T : TFVec
        Temperature [C]
    T : TFVec
        Salinity [ppm]
    model : str
        Model used to calculate solubility. Defaults to 'rpl'.

    Returns
    -------
    TFVec
        Solubility of petroleum gas in brine.
    """
    v = BRINE.call(model, "gas_solubility", P, T, Sa)
    if isinstance(P, np.ndarray):
        return np.array(v)
    return v
