# Instructions to run Particle-In-Cell (PIC) simulations

This document explains how to perform PIC simulations with a 1D Python code and with the Smilei code.

First, you need to prepare your laptop for the activity. You should have a Linux machine (virtual if you use Windows) or MacOS with a running Python interpreter and the PIC code [Smilei](https://smileipic.github.io/Smilei/) installed.
This [guide](./smilei_guide.md)

- Prepare an input file ``input.py`` with this content.

```python
########################################
## ---------------------------------- ##
## INPUT FILE FOR THE PIC-CODE SMILEI ##
##       1D PIC LPP SIMULATIONS       ##
## ---------------------------------- ##
########################################
import math
import numpy as np
from math import pi, sqrt, sin, tan, log10, exp, ceil, log, atan
from scipy.constants import pi, c, m_e, h, hbar

## ---------------------------------- ##
##             PARAMETERS             ##
## ---------------------------------- ##

# UNITS
lambda_SI = 0.8e-6
omega_SI = 2.*pi*c/lambda_SI
micron = um = 1.e-6  * omega_SI/c
fs = 1.e-15 * omega_SI
pemr = 1836.15267343

# DISCRETIZATION
dx = um / 20
Lx = 60 * um
nx = int( Lx / dx )
Tsim = 300 * fs
cfl = 0.98
dt = cfl * dx
nt = int(Tsim / dt)

# PLASMA PARAMETERS
Zp = 1.
Ap = 1.
ne0 = 1e-2
nppc = 2
me = 1
mi = pemr * Ap

# LASER PARAMETERS
a0 = 1
fwhm = 30 * fs
laser_center = 2 * fwhm

# DIAGNOSTICS
every_fs = int(fs/dt)
every_out = 20*every_fs
shift = every_out
## ---------------------------------- ##
##         SIMULATION SETUP           ##
## ---------------------------------- ##

Main(
    geometry = "1Dcartesian",
    cell_length = [dx],
    grid_length = [Lx],
    timestep = dt,
    simulation_time = Tsim,
    interpolation_order = 2,
    EM_boundary_conditions = [ ['silver-muller'] ],
    maxwell_solver = 'Yee',
    random_seed = smilei_mpi_rank,
    print_every = int(nt/100.0),
    reference_angular_frequency_SI = omega_SI,
    solve_poisson = False,
    number_of_patches = [2],
    patch_arrangement = 'hilbertian',
)

## ---------------------------------- ##
##                LASER               ##
## ---------------------------------- ##

LaserPlanar1D(
  box_side = 'xmin',
  a0 = a0,
  omega = 1,
  polarization_phi = 0, 
  ellipticity = 0,
  time_envelope = tgaussian(fwhm = fwhm, center = laser_center)
)

## ---------------------------------- ##
##               PLASMA               ##
## ---------------------------------- ##

# ELECTRONS
Species(
  name = "ELE",
  position_initialization = "regular",
  regular_number = [nppc],
  momentum_initialization = "cold",
  particles_per_cell = nppc,
  mass = me,
  number_density = trapezoidal(ne0, xvacuum = 0.4 * Lx, xplateau = 0.2 * Lx, xslope1 = 0, xslope2 = 0), 
  charge = -1.,
  boundary_conditions = [["reflective", "reflective"]],
  pusher = "boris",
)

# IONS
Species(
  name = "ION",
  position_initialization = "regular",
  regular_number = [nppc],
  momentum_initialization = "cold",
  particles_per_cell = nppc,
  mass = mi,
  number_density = trapezoidal(ne0 / Zp, xvacuum = 0.4 * Lx, xplateau = 0.2 * Lx, xslope1 = 0, xslope2 = 0), 
  charge = Zp,
  boundary_conditions = [["reflective", "reflective"]],
  pusher = "boris",
)

## ---------------------------------- ##
##            DIAGNOSTICS             ##
## ---------------------------------- ##

# SCALARS
DiagScalar(every=1.0)

DiagTrackParticles(
    species = "ELE",
    every = [shift, every_out], # [start, period]
    attributes = ["x", "w", "px", "py", "pz"]
)

# TRACKPARTICLES
DiagTrackParticles(
    species = "ION",
    every = [shift, every_out],
    attributes = ["x", "w", "px", "py", "pz"]
)

# FIELDS
DiagFields(
    every = [shift, every_out],
    time_average = 1,
    fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jx", "Jy", "Jz", "Rho_ELE", "Rho_ION", "Rho"],
)
```
- Create a folder named smilei (``mkdir smilei``) inside the directory ``Lecture2``. 

- Copy there the input file.

- Run a test simulation with this input file following the instructions at the end of [Smilei Guide](../smilei_guide.md):
  
  type in a terminal
  ```
  export OMP_NUM_THREADS=2
  export OMP_SCHEDULE=dynamic
  ```
  and run the simulation
  ```
  mpirun -np 2 <path_to_smilei_folder>/smilei input.py
  ```

- In the directory ``Lecture2`` run the script [plot.py](./plot.py) by typing in the terminal ``python plot.py``. This script prints plots comparing the results of this code and the results of the 1D Python PIC code in the ``Data`` folder. Check to have installed all the required Python modules.

- Now you can visualise some plots in the folder ``Output``.