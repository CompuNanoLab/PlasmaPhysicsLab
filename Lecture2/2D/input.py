########################################
## ---------------------------------- ##
## INPUT FILE FOR THE PIC-CODE SMILEI ##
##    2D LASER-PLASMA INTERACTION     ##
## ---------------------------------- ##
########################################
import math
import numpy as np
from math import pi, sqrt, sin, tan, log10, exp, ceil, log, atan
from scipy.constants import pi, c, m_e, h, hbar

## ---------------------------------- ##
##             PARAMETERS             ##
## ---------------------------------- ##

# GENERAL PARAMETERS
lambda_SI           = 0.8e-6
omega_SI            = 2.*pi*c/lambda_SI
um                  = 1.e-6  * omega_SI/c
fs                  = 1.e-15 * omega_SI
pemr                = 1836.15267343
mc2                 = 0.510998950e6
mev                 = 1/0.510998950

# BOX PARAMETERS
target_thickness    = 20*um
targetstart         = 20*um
aftertarget         = 20*um
Lx                  = targetstart+target_thickness+aftertarget
Ly                  = 40*um
resx                = 20
resy                = 20
dx                  = um/resx
dy                  = um/resy
nx                  = Lx/dx
ny                  = Ly/dy
npatch_x            = 1
npatch_y            = 4
Lx                  = ceil(nx/npatch_x)*npatch_x*dx
Ly                  = ceil(ny/npatch_y)*npatch_y*dy

#TIME
Tsim                = Lx + 60*fs
cfl                 = 0.98
dt                  = cfl/sqrt(1./dx**2+1./dy**2)
Nsteps              = int(Tsim/dt)

# LASER PARAMETERS
a0                  = 40
waist               = 2.8*um/sqrt(2*log(2))
incidence_angle     = 0
pulse_duration      = 23*fs
offset              = targetstart-Ly/2*sin(incidence_angle)
def pulse_temporal_profile(t):
    if 0 < t < 2*pulse_duration:
        return sin(pi/2. * t/pulse_duration)
    else : return 0.

# PLASMA PARAMETERS
Temp                = 10/mc2
Z_C                 = 6
A_C                 = 12
species_boundary_conditions=[['remove'],['periodic'],]
n_target            = 5
nppc_ion            = 1
nppc_ele            = 2

# DIAGNOSTIC PARAMETERS
every_fs            = int(fs/dt)
timescalar          = [0,  int(0.5*every_fs)]
timefield           = [10*every_fs,int(Tsim/dt),30*every_fs]
timetrack           = [10*every_fs,int(Tsim/dt),30*every_fs]

## ---------------------------------- ##
##         SIMULATION SETUP           ##
## ---------------------------------- ##

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2 ,
    cell_length = [dx,dy],
    grid_length  = [Lx,Ly],
    number_of_patches = [npatch_x,npatch_y],
    timestep = dt,
    simulation_time = Tsim,
    reference_angular_frequency_SI = omega_SI,
    EM_boundary_conditions = [['silver-muller'],['periodic']],
    solve_poisson= False,
    patch_arrangement = 'hilbertian',
    maxwell_solver = 'Yee',
    random_seed = smilei_mpi_rank,
    print_every = 40,
)

## ---------------------------------- ##
##                LASER               ##
## ---------------------------------- ##

LaserGaussian2D(
    box_side            = 'xmin',
    a0                  = a0 ,
    waist               = waist,
    focus               = [targetstart, Ly / 2],
    ellipticity         = 0.,
    incidence_angle     = incidence_angle,
    time_envelope       = pulse_temporal_profile
) 

## ---------------------------------- ##
##               PLASMA               ##
## ---------------------------------- ##

# IONS
Species(
    name = "ion",
    position_initialization = 'random',
    particles_per_cell = nppc_ion,
    momentum_initialization = "cold",
    number_density = polygonal(xpoints=[targetstart, targetstart+target_thickness],xvalues=[n_target/Z_C, n_target/Z_C]),
    mass = A_C*pemr,
    charge = Z_C,
    boundary_conditions = species_boundary_conditions,
    pusher = "boris"
)

# ELECTRONS
Species(
    name = "ele",
    position_initialization = "random",
    particles_per_cell = nppc_ele,
    momentum_initialization = "maxwell-juettner",
    number_density =polygonal(xpoints=[targetstart, targetstart+target_thickness],xvalues=[n_target, n_target]), 
    temperature = [Temp],
    mean_velocity = [0.],
    charge = -1.,
    mass = 1.,
    boundary_conditions = species_boundary_conditions,
    pusher = "boris",
)

## ---------------------------------- ##
##            DIAGNOSTICS             ##
## ---------------------------------- ##

# SCALARS
DiagScalar(
    every = timescalar
)

# FIELDS
DiagFields(
    every = timefield,
    flush_every=1000,
    fields = ["Rho_ele","Rho_ion","Bx","By","Bz","Ex","Ey","Ez"]
)

# BINNING - PARTICLES
#ions
DiagParticleBinning(
    deposited_quantity = "weight",
    every = timescalar,
    time_average = 1,
    species = ["ion"],
    axes = [ ["ekin",    0.,	200.*mev,   2000] ]
)

#ele
DiagParticleBinning(
    deposited_quantity = "weight",
    every = timescalar,
    time_average = 1,
    species = ["ele"],
    axes = [ ["ekin",    0.,	300.*mev,   3000] ]
)
