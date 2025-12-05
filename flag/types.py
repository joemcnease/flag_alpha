"""
Custom types for flag package.
"""
from typing import TypeVar

import numpy as np


TFVec = TypeVar("FVec", float, list[float], np.ndarray)
