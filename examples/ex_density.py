import numpy as np
import matplotlib.pyplot as plt

import flag


Ts    = [120, 100, 80, 60, 40, 20]  # C
P     = np.linspace(5, 60, 100)     # MPa
Sa    = np.ones_like(P)*200000      # ppm
NaCl  = np.ones_like(P)*100         # weight %
KCl   = np.zeros_like(P)            # weight %
CaCl2 = np.zeros_like(P)            # weight %

fig, axs = plt.subplots(1, 2)

for i in range(len(Ts)):
    T = np.ones_like(P)*Ts[i]

    rho_h2o   = flag.h2o.density(P, T)
    rho_brine = flag.brine.density(P, T, Sa, NaCl, KCl, CaCl2)

    axs[0].plot(P, rho_h2o, label=f"T = {Ts[i]}")
    axs[1].plot(P, rho_brine)

axs[0].set_ylabel("Density [g/cm^3]")
titles = ["H2O", "Brine"]
for i in range(2):
    axs[i].set_xlabel("Pressure [MPa]")
    axs[i].set_xlim(P.min(), P.max())
    axs[i].set_title(titles[i])
    axs[i].grid()

fig.legend(loc='upper right', bbox_to_anchor=(1.0, 0.89))

plt.show()
