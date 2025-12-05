# flag

![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![Documentation](https://img.shields.io/badge/docs-in_docs%2F-brightgreen)
![Build](https://img.shields.io/badge/build-passing-brightgreen)

flag is a Python package for computing **fluid properties** such as **density, velocity, bulk modulus, viscosity,** and more. These properties are useful for **oil exploration** and **CCUS (Carbon Capture, Utilization, and Storage) monitoring projects**.  

The package currently includes modules for **H2O**, **Brine**, **CO2**, **Gas**, **Oil**, and more. See
the full documentation in docs/ for detailed overview of all available modules.

Additional modules are under development to expand functionality.

---

## Features

- Computes key fluid properties for engineering and geoscience applications  
- Supports multiple fluids including water, brine, CO2, gas, oil, and more
- Works on Windows, macOS, and Linux  
- Supports Python 3.9 and above  

---

## Installation

### Using Pre-built Wheels

We provide pre-built wheels for different operating systems and Python versions. If you have the wheels in a folder, you can install using `pip`:

```bash
pip install --no-index --find-links=/path/to/flagpy_wheels flag
```

Replace /path/to/flagpy_wheels with the path to your wheel folder.

Pip will automatically select the correct wheel for your Python version and OS.

## Quick Example

```python
import numpy as np
import matplotlib.pyplot as plt

import flag


Ts    = [120, 100, 80, 60, 40, 20]  # C
P     = np.linspace(5, 60, 100)     # MPa
Sa    = np.ones_like(P)*200000      # ppm
NaCl  = np.ones_like(P)*100         # weight %
KCl   = np.zeros_like(P)            # weight %
CaCl2 = np.zeros_like(P)            # weight %

fig, axs = plt.subplots(1, 3)

for i in range(len(Ts)):
    T = np.ones_like(P)*Ts[i]

    rho_h2o   = flag.h2o.density(P, T)
    rho_brine = flag.brine.density(P, T, Sa, NaCl, KCl, CaCl2)
    rho_co2   = flag.co2.density(P, T)

    axs[0].plot(P, rho_h2o, label=f"T = {Ts[i]}")
    axs[1].plot(P, rho_brine)
    axs[2].plot(P, rho_co2)

axs[0].set_ylabel("Density [g/cm^3]")
titles = ["H2O", "Brine", "CO2"]
for i in range(3):
    axs[i].set_xlabel("Pressure [MPa]")
    axs[i].set_xlim(P.min(), P.max())
    axs[i].set_title(titles[i])
    axs[i].grid()

fig.legend(loc='upper right', bbox_to_anchor=(1.0, 0.89))

plt.show()
```

You can also refer to the Getting Started guide for more examples and usage instructions.

Supported Platforms
    Windows (Python >= 3.9)
    macOS (Python >= 3.9)
    Linux (Python >= 3.9)

If you require support for other platforms or Python versions, please contact us.

## Documentation

Full documentation is available in the docs/ folder.


## Contact

For questions, issues, or suggestions, please reach out:

    Email: Yingcai Zheng – yzheng24@central.uh.edu
    Email: Min Sun       – msun@central.uh.edu
    Email: Joe McNease   – jmcnease@uh.edu
