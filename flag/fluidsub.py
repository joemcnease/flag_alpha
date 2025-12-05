"""
Some helpful tools for fluid substitution.

See Gassman (1951) and Smith, Sondergeld, and Rai (2003)
for more information on theory and practice, respectively.
"""
from .types import TFVec

import numpy as np


def fluidsub(Ks: TFVec, K0: TFVec, Kfl: TFVec, phi: TFVec):
    """
    Isotropic fluid subsitution using low-frequency
    Gassman (1951) theory.

    Parameters
    ----------
    Ks : TFVec
        Bulk modulus of rock frame [GPa]
    K0 : TFVec
        Bulk modulus of mineral matrix [GPa]
    Kfl : TFVec
        Bulk modulus of pore fluid [GPa]
    phi : TFVec
        Porosity

    Returns
    -------
    TFVec
        Saturated bulk modulus
    """
    num = (K0**2)*(1.0 - Ks/K0)**2
    den = (K0**2/Kfl)*phi + K0*(1.0 - phi) - (K0*Ks)**2

    return Ks + num/den


def reuss_average(K: TFVec, S: TFVec):
    """
    Reuss (isostress) average.

    Parameters
    ----------
    K : TFVec
        Quantity to average.
    S : TFVec
        Averaging weights for each quantity.

    Returns
    -------
    TFVec
        Reuss average of K vector
    """
    return 1.0/np.sum(S/K)


def voigt_average(K: TFVec, S:TFVec):
    """
    Voigt (isostrain) average.

    Parameters
    ----------
    K : TFVec
        Quantity to average.
    S : TFVec
        Averaging weights for each quantity.

    Returns
    -------
    TFVec
        Voigt average of K vector
    """
    return np.sum(K*S)


def vrh_average(K: TFVec, S: TFVec):
    """
    Voigt-Reuss-Hill average.

    Parameters
    ----------
    K : TFVec
        Quantity to average.
    S : TFVec
        Averaging weights for each quantity.

    Returns
    -------
    TFVec
        Voigt-Reuss-Hill average of K vector
    """
    return 0.5*(reuss_average(K, S) + voigt_average(K, S))
