import numpy as np
import matplotlib.pyplot as plt

import flag


P = 30
T = 15
V = flag.h2o.velocity(P, T)
print(V)

P = [30, 30]
T = [15, 15]
V = flag.h2o.velocity(P, T)
print(V)

P = np.array([30, 30])
T = np.array([15, 15])
V = flag.h2o.velocity(P, T)
print(V)
