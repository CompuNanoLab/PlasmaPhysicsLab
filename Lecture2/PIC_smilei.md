# Instructions to run a PIC simulation with Smilei

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
Lx = 50 * um 
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
shift = int(4*fwhm/dt)
every_fs = int(fs/dt)
every_out = 20*every_fs

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
  position_initialization = "random",
  momentum_initialization = "cold",
  particles_per_cell = nppc,
  mass = me,
  number_density = trapezoidal(n0, xvacuum = 0.4 * Lx, xplateau = 0.6 * Lx, xslope1 = 0, xslope2 = 0), 
  charge = -1.,
  boundary_conditions = [["periodic", "periodic"]],
  pusher = "boris",
)

# IONS
Species(
  name = "ION",
  position_initialization = "random",
  momentum_initialization = "cold",
  particles_per_cell = nppc,
  mass = mi,
  number_density = trapezoidal(n0, xvacuum = 0.4 * Lx, xplateau = 0.6 * Lx, xslope1 = 0, xslope2 = 0), 
  charge = Z,
  boundary_conditions = [["periodic", "periodic"]],
  pusher = "boris",
)

## ---------------------------------- ##
##            DIAGNOSTICS             ##
## ---------------------------------- ##

# SCALARS
DiagScalar(every=1.0)

DiagTrackParticles(
    species = "ele",
    every = [shift, every_out], # [start, period]
    attributes = ["x","w", "px", "py", "pz"]
)

# TRACKPARTICLES
DiagTrackParticles(
    species = "ion",
    every = [shift, every_out],
    attributes = ["x","w", "px", "py", "pz"]
)

# FIELDS
DiagFields(
    every = [shift, every_out],
    time_average = 1,
    fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jx_ELE", "Jy_ELE", "Jz_ELE", "Jx_ION", "Jy_ION", "Jz_ION", "Rho_ELE", "Rho_ION", "Rho"],
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
  mpirun -np 4 <path_to_smilei_folder/smilei input.py
  ```

- In the directory ``Lecture2`` run the script ``plot.py`` by typing in the terminal ``python plot.py``.

- Now you can visualise some plots in the folder ``output``.
