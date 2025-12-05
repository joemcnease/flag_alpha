Getting Started
===============

This section will help you get started with the `flag` package.

Installation
------------

It is generally recommended to install Python packages inside a virtual
environment so that dependencies for different projects do not interfere
with one another. If you are new to Python, here are two common ways to
create an isolated environment.

Creating a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tabs::

   .. tab:: Using Conda

      .. code-block:: bash

          conda create -n flag-env python=3.11
          conda activate flag-env

   .. tab:: Using venv

      .. code-block:: bash

          python3 -m venv env
          source env/bin/activate   # Linux/Mac
          env\Scripts\activate      # Windows


You can install the package using pip:

.. code-block:: bash

    pip install --no-index --find-links=/path/to/flagpy_wheels flag

where /path/to/flagpy_wheels is the path to the directory which contains
the pre-built wheels for different python versions and operating systems.
We currently have pre-built wheels for Windows, macOS, and Linux for
Python versions >= 3.9. If you have other system or Python version requirements,
please see the :ref:`contact` page.

Quick Example
-------------

Here’s a simple example using the package:

.. code-block:: python

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
    for i in range(3):
        axs[i].set_xlabel("Pressure [MPa]")
        axs[i].set_xlim(P.min(), P.max())
        axs[i].grid()
    
    fig.legend(loc='upper right', bbox_to_anchor=(1.0, 0.89))
    
    plt.show()

.. image:: _static/density_example.png
   :alt: Quick example figure
   :align: center
   :width: 800px

Next Steps
----------

- Read the Modules section for detailed function/class usage.
- Explore examples in the `examples/` folder.
