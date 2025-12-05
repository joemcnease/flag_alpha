import numpy as np
import matplotlib.pyplot as plt

import flag


def main():
    # For a single input
    Sa    = 200000
    NaCl  = 100.0
    KCl   = 0.0
    CaCl2 = 0.0
    T     = 20.0
    P     = 10.0
    Vb = flag.brine.velocity(P, T, Sa, NaCl, KCl, CaCl2)

    # For lists
    # For a single input
    Sa    = [200000]
    NaCl  = [100.0]
    KCl   = [0.0]
    CaCl2 = [0.0]
    T     = [20.0]
    P     = [10.0]
    Vb = flag.brine.velocity(P, T, Sa, NaCl, KCl, CaCl2)

    # For numpy arrays (and a range of pressures and temperatures)
    fig, axs = plt.subplots(1, 1)
    for i in range(2, 10):
        Sa    = np.ones(100)*200000
        NaCl  = np.ones(100)*100.0
        KCl   = np.ones(100)*0.0
        CaCl2 = np.ones(100)*0.0
        T     = np.ones(100)*(i*30)
        P     = np.linspace(10, 100, 100)
        Vb = flag.brine.velocity(P, T, Sa, NaCl, KCl, CaCl2)

        axs.plot(P, Vb, label=f"T={np.round(T[0])}")

    axs.set_xlabel("Pressure [MPa]")
    axs.set_ylabel("Velocity Brine [km/s]")
    axs.grid()
    axs.legend()
    plt.show()

    # P-T surface
    P = np.linspace(5, 100)
    T = np.linspace(5, 100)
    Sa    = np.ones_like(T)*200000
    NaCl  = np.ones_like(T)*100.0
    KCl   = np.ones_like(T)*0.0
    CaCl2 = np.ones_like(T)*0.0
    Vb = np.zeros([P.shape[0], T.shape[0]])
    for i in range(len(P)):
        Ptmp = np.ones_like(T)*P[i]
        Vb[i, :] = flag.brine.velocity(Ptmp, T, Sa, NaCl, KCl, CaCl2)

    fig, axs = plt.subplots(1, 1)
    axs.contourf(Vb, origin="lower")
    axs.set_xlabel("Temperature [C]")
    axs.set_ylabel("Pressure [MPa]")
    plt.show()

    # Velocity, density, and bulk modulus
    T = np.linspace(5, 100)
    Sa    = np.ones_like(T)*0
    NaCl  = np.ones_like(T)*100.0
    KCl   = np.ones_like(T)*0.0
    CaCl2 = np.ones_like(T)*0.0

    fig, axs = plt.subplots(1, 3)
    ylbls = ["Velocity [km/s]", "Density [g/cm^3]", "Bulk Modulus [GPa]"]

    Ps = [0, 40, 80, 120]
    for i in range(len(Ps)):
        P = np.ones_like(T)*Ps[i]
        Vb = flag.brine.velocity(P, T, Sa, NaCl, KCl, CaCl2)
        rhob = flag.brine.density(P, T, Sa, NaCl, KCl, CaCl2)
        Kb = flag.brine.bulk_modulus(P, T, Sa, NaCl, KCl, CaCl2)
        vals = [Vb, rhob, Kb]
        for j in range(len(vals)):
            axs[j].plot(T, vals[j], label=f"P [MPa] = {Ps[i]:d}")
            axs[j].set_xlabel("Temperature [C]")
            axs[j].set_ylabel(ylbls[j])

    [axs[i].grid() for i in range(len(axs))]
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys())
    plt.show()


if __name__ == "__main__":
    main()
