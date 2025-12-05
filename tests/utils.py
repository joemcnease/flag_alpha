from typing import TypeVar

import numpy as np


TFVec = TypeVar("TFVec", float, list[float], np.ndarray)


def close_to(val1: float, val2: float, tol: float = 1e-4) -> bool:
    """
    Compare equivalence of floating point values with a given tolerance.
    """
    return abs(val1-val2) < tol


def close_to_array(lst1: list[float], lst2: list[float], tol: float = 1e-4) -> bool:
    """
    Compare equivalence of floating point arrays (list, numpy array, etc.)
    with a given tolerance.
    """
    assert len(lst1) == len(lst2)
    res = sum([close_to(lst1[i], lst2[i], tol) for i in range(len(lst1))])
    return res == len(lst1)


def to_API(rho0: TFVec) -> TFVec:
    """
    Compute API from density at STP.
    """
    return 141.5/rho0 - 131.5


def from_API(api: TFVec) -> TFVec:
    """
    Compute density at STP from API.
    """
    return 1./((api + 131.5)/141.5)
