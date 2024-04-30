# Instruction to run a PIC simulation with Smilei

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
ZC = 6.
AC = 6.
ne0 = 1
nppc = 2
me = 1
mi = pemr * AC

p0x=0.
p0y=0.
p0z=0.
dpx=0.
dpy=0.
dpz=0.
g0 = sqrt(1.+p0x**2+p0y**2+p0z**2)
v0ex = p0x/g0/me
v0ey = p0y/g0/me
v0ez = p0z/g0/me
Tx = dpx**2/me #mass*dpx**2
Ty = dpy**2/me #mass*dpy**2
Tz = dpz**2/me #mass*dpz**2
Te = [dpx**2,dpy**2,dpz**2]
#Te = [Tx, Ty, Tz]
#Ti = [Tx, Ty, Tz]

g0 = sqrt(1.+(p0x**2+p0y**2+p0z**2)/mi**2)
v0ix = p0x/g0/mi
v0iy = p0y/g0/mi
v0iz = p0z/g0/mi
Ti = [dpx**2/mi,dpy**2/mi,dpz**2/mi]

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
    interpolation_order = 2,
    cell_length = [dx],
    #number_of_cells = [nx], 
    grid_length = [Lx],
    number_of_patches = [2],
    timestep = dt,
    simulation_time = T_sim,
    EM_boundary_conditions = [ ['silver-muller'] ],
    random_seed = smilei_mpi_rank,
    print_every = int(nt/100.0),
    reference_angular_frequency_SI = omega_SI,
    solve_poisson = False
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
  momentum_initialization = "rectangular",
  mean_velocity = [v0ex, v0ey, v0ez],
  temperature = Te,
  particles_per_cell = nppc,
  mass = me,
  number_density = trapezoidal(n0, xvacuum = 3 * fwhm, xplateau = Lx - 2 * fwhm, xslope1 = 0, xslope2 = 0), 
  charge = -1.,
  boundary_conditions = [["reflective", "reflective"]],
  time_frozen = 0.0,
  is_test = False,
  pusher = "boris",
)

# IONS
Species(
    name = "ION",
    position_initialization = "regular",
    regular_number = [nppc],
    momentum_initialization = "rectangular",
    mean_velocity = [v0ix, v0iy, v0iz],
    temperature = Ti,
    particles_per_cell = nppc,
    mass = mi,
    number_density = trapezoidal(n0 / Z, xvacuum = 3 * fwhm, xplateau = Lx - 2 * fwhm, xslope1 = 0, xslope2 = 0), 
    charge = Z,
    boundary_conditions = [["reflective", "reflective"]],
    time_frozen = 0.0,
    is_test = False,
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
    fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jx_ele", "Jy_ele", "Jz_ele","Jx_ion", "Jy_ion", "Jz_ion", "Rho_ele", "Rho_ion", "Rho"],
)

```
- Create a folder named ``smilei`` inside the directory ``Lecture2``. 

- Copy there the input file.

- Run a test simulation with this input file following the instructions at the end of [Smilei Guide](../smilei_guide.md).

- In the directory ``Lecture2`` run the script ``plot.py`` by typing in a terminal ``python plot.py``.

- Now you can visualise some plots in the folder ``output``.
