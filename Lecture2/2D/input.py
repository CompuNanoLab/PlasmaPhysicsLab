########################################
## ---------------------------------- ##
## INPUT FILE FOR THE PIC-CODE SMILEI ##
##    2D LASER-PLASMA INTERACTION     ##
## ---------------------------------- ##
########################################

import math
import numpy as np
from math import pi, sqrt, sin, log, ceil
from scipy.constants import micron, c, pi, m_e

## ---------------------------------- ##
##             PARAMETERS             ##
## ---------------------------------- ##

# GENERAL PARAMETERS
lambda_SI = 0.8e-6  # Wavelength of the laser in meters
omega_SI = 2. * pi * c / lambda_SI  #Laser frequency
um = 1.e-6 * omega_SI / c  # Conversion factor from meters to micrometers
fs = 1.e-15 * omega_SI  # Conversion factor from seconds to femtoseconds
pemr = 1836.15267343  # Electron-proton mass ratio
mc2 = 0.510998950e6  # Electron rest mass energy in eV
mev = 1 / 0.510998950  # Conversion factor from eV to MeV

# BOX PARAMETERS
target_thickness = 20 * um  # Thickness of the target
targetstart = 20 * um  # Start position of the target
aftertarget = 20 * um  # Distance after the target
Lx = targetstart + target_thickness + aftertarget  # Total length of the box along x-axis
Ly = 40 * um  # Total length of the box along y-axis
resx = 20  # Points per micron in x-direction
resy = 20  # Points per micron in x-direction
dx = um / resx  # Cell size in x-direction
dy = um / resy  # Cell size in y-direction
nx = Lx / dx  # Number of cells in x-direction
ny = Ly / dy  # Number of cells in y-direction
npatch_x = 1  # Number of patches in x-direction for workload distribution
npatch_y = 4  # Number of patches in y-direction
Lx = ceil(nx / npatch_x) * npatch_x * dx  # Adjusted length of the box along x-axis
Ly = ceil(ny / npatch_y) * npatch_y * dy  # Adjusted length of the box along y-axis

# TIME
Tsim = 2*Lx # Total simulation time
cfl = 0.98  # CFL coefficient
dt = cfl / sqrt(1. / dx**2 + 1. / dy**2)  # Time step size
Nsteps = int(Tsim / dt)  # Number of time steps

# LASER PARAMETERS
a0 = 40  # Laser amplitude
waist = 2.8 * um / sqrt(2 * log(2))  # Waist of the laser beam
incidence_angle = 0  # Incidence angle of the laser
pulse_duration = 23 * fs  # Duration of the laser pulse
offset = targetstart - Ly / 2 * sin(incidence_angle)  # Offset of the laser pulse
# Temporal profile of the laser pulse
def pulse_temporal_profile(t):
    if 0 < t < 2 * pulse_duration:
        return sin(pi / 2. * t / pulse_duration)
    else:
        return 0.

# PLASMA PARAMETERS
Temp = 10 / mc2  # Temperature of the plasma in MeV
Z_C = 6  # Charge number of the ions (Carbon)
A_C = 12  # Atomic mass number of the ions (Carbon)
species_boundary_conditions = [['remove'], ['periodic']]  # Boundary conditions for species
n_target = 5  # Number density of the target in units of critical density
nppc_ion = 1  # Number of ions per cell
nppc_ele = 2  # Number of electrons per cell

# DIAGNOSTIC PARAMETERS
every_fs = int(fs / dt)  # Frequency of diagnostics in terms of femtoseconds
timescalar = [0, int(0.5 * every_fs)]  # Time steps for scalar diagnostics
timefield = [10 * every_fs, int(Tsim / dt), 30 * every_fs]  # Time steps for field diagnostics
timetrack = [10 * every_fs, int(Tsim / dt), 30 * every_fs]  # Time steps for particle tracking diagnostics

## ---------------------------------- ##
##         SIMULATION SETUP           ##
## ---------------------------------- ##

Main(
    geometry="2Dcartesian",
    interpolation_order=2,
    cell_length=[dx, dy],
    grid_length=[Lx, Ly],
    number_of_patches=[npatch_x, npatch_y],
    timestep=dt,
    simulation_time=Tsim,
    reference_angular_frequency_SI=omega_SI,
    EM_boundary_conditions=[['silver-muller'], ['periodic']],
    solve_poisson=False,
    patch_arrangement='hilbertian',
    maxwell_solver='Yee',
    random_seed=smilei_mpi_rank,
    print_every=40,
)

## ---------------------------------- ##
##                LASER               ##
## ---------------------------------- ##

LaserGaussian2D(
    box_side='xmin',
    a0=a0,
    waist=waist,
    focus=[targetstart, Ly / 2],
    ellipticity=0.,
    incidence_angle=incidence_angle,
    time_envelope=pulse_temporal_profile
)

## ---------------------------------- ##
##               PLASMA               ##
## ---------------------------------- ##

# IONS
Species(
    name="ion",
    position_initialization='random',
    particles_per_cell=nppc_ion,
    momentum_initialization="cold",
    number_density=polygonal(xpoints=[targetstart, targetstart + target_thickness],
                             xvalues=[n_target / Z_C, n_target / Z_C]),
    mass=A_C * pemr,
    charge=Z_C,
    boundary_conditions=species_boundary_conditions,
    pusher="boris"
)

# ELECTRONS
Species(
    name="ele",
    position_initialization="random",
    particles_per_cell=nppc_ele,
    momentum_initialization="maxwell-juettner",
    number_density=polygonal(xpoints=[targetstart, targetstart + target_thickness], xvalues=[n_target, n_target]),
    temperature=[Temp],
    mean_velocity=[0.],
    charge=-1.,
    mass=1.,
    boundary_conditions=species_boundary_conditions,
    pusher="boris",
)

## ---------------------------------- ##
##            DIAGNOSTICS             ##
## ---------------------------------- ##

# SCALARS
DiagScalar(
    every=timescalar
)

# FIELDS
DiagFields(
    every=timefield,
    flush_every=1000,
    fields=["Rho_ele", "Rho_ion", "Bx", "By", "Bz", "Ex", "Ey", "Ez"]
)

# BINNING - PARTICLES
# Ions
DiagParticleBinning(
    deposited_quantity="weight",
    every=timescalar,
    time_average=1,
    species=["ion"],
    axes=[["ekin", 0., 200. * mev, 2000]]
)

# Electrons
DiagParticleBinning(
    deposited_quantity="weight",
    every=timescalar,
    time_average=1,
    species=["ele"],
    axes=[["ekin", 0., 300. * mev, 3000]]
)
