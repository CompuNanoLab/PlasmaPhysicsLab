# Instructions to run Particle-In-Cell (PIC) simulations

This document explains you how to perform PIC simulations with a 1D Python code and with the Smilei code as done during the lecture on 02/05/2024.

## Before starting

First, you need to prepare your laptop for the activity. You should have a Linux machine (virtual if you use Windows) or MacOS with a running Python interpreter and the PIC code [Smilei](https://smileipic.github.io/Smilei/) installed.

This [guide](./smilei_guide.md) helps you to prepare the machine. 

Here is a list of basic useful commands that you can use in a terminal:

```
    cat --- for creating and displaying short files
    chmod --- change permissions
    cd --- change directory
    clear --- clear screen
    cp --- for copying files
    diff --- compares files
    echo --- echo argument
    grep --- search file
    head --- display first part of file
    history --- show history of previous commands
    ls --- see what files you have
    less --- use to read files (q to exit)
    man --- view manual pages for Unix commands
    more --- use to read files
    mkdir --- create directory
    mv --- for moving and renaming files
    pwd --- find out what directory you are in
    rm --- remove a file
    rmdir --- remove directory
    setenv --- set an environment variable
    sort --- sort file
    tail --- display last part of file
    tar --- create an archive, add or extract files
    wc --- count characters, words, lines
    who --- tells you who's logged on

```

To edit files you can use a command-line editor like Nano. To open a file using Nano, open your terminal and type the following command:

```bash
nano text.txt
```
Then, you can move with arrow keys and edit the file as you want. Use Ctrl + O to save changes to the file and Ctrl + X to exit.

## Running simulations

Let's start! open a terminal in your machine and do the following steps.
Download the GitHub repository containing all this material:

```bash
git clone https://github.com/CompuNanoLab/PlasmaPhysicsLab.git
```

Go to the right folder:

```bash
cd PlasmaPhysicsLab/Lecture2
```

Check that you have installed what you need to run Jupyter Notebooks (you may need to use pip3 instead of pip):

```bash
pip install jupyter
```

Open the notebook with the 1D Python PIC code:

```bash
jupyter-notebook 1DPIC.ipynb
```

a browser page should start where you can edit and run all the cells in the notebook. The output of the code is saved in the ```Data``` folder.

Now open a new terminal and move again to the ```Lecture2``` folder:

```bash
cd PlasmaPhysicsLab/Lecture2
```

Generate a folder in which you will perform the Smilei simulation:

```bash
mkdir smilei
cd smilei
```
Prepare an input file ``input.py`` for the simulation:

```bash
nano input.py
```
where you can copy the following lines:

```python
########################################
## ---------------------------------- ##
## INPUT FILE FOR THE PIC-CODE SMILEI ##
##       1D PIC LPP SIMULATIONS       ##
## ---------------------------------- ##
########################################

# Importing necessary libraries
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
pemr = 1836.15267343  # Electron-to-Proton mass ratio

# DISCRETIZATION
dx = um / 20  # Grid spacing
Lx = 60 * um  # Total length of the domain
nx = int( Lx / dx )  # Number of grid points
Tsim = 300 * fs # Total simulation time
cfl = 0.98  # CFL number
dt = cfl * dx  # Time step size
nt = int(Tsim / dt)  # Number of time steps

# PLASMA PARAMETERS
Zp = 1.  # Ion charge
Ap = 1.  # Ion atomic number
ne0 = 1e-2  # Electron density
nppc = 2  # Number of particles per cell
me = 1  # Electron mass
mi = pemr * Ap  # Ion mass

# LASER PARAMETERS
a0 = 1  # Peak value of the normalized vector potential
fwhm = 30 * fs  # Full width half maximum in intensity
laser_center = 2 * fwhm  # Temporal centering of the laser

# DIAGNOSTICS
every_fs = int(fs/dt)
every_out = 20*every_fs  # Output interval
shift = every_out  # When to start printing diagnostics

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
    number_of_patches = [2],  # Patches for workload distribution
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
  time_envelope = tgaussian(fwhm = fwhm, center = laser_center)  # Gaussian time envelope
)

## ---------------------------------- ##
##               PLASMA               ##
## ---------------------------------- ##

# ELECTRONS
Species(
  name = "ELE",
  position_initialization = "regular", # Particles are distributed equally spaced in the cell
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
DiagScalar(every=1.0)  # Scalar diagnostics

# TRACKPARTICLES
DiagTrackParticles(
    species = "ELE",
    every = [shift, every_out],  # Output interval for electron species
    attributes = ["x", "w", "px", "py", "pz"]  # Attributes to track for electrons
)

# TRACKPARTICLES
DiagTrackParticles(
    species = "ION",
    every = [shift, every_out],  # Output interval for ion species
    attributes = ["x", "w", "px", "py", "pz"]  # Attributes to track for ions
)

# FIELDS
DiagFields(
    every = [shift, every_out],  # Output interval for fields
    time_average = 1,  # Time averaging
    fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jx", "Jy", "Jz", "Rho_ELE", "Rho_ION", "Rho"],  # Fields to save
)
```

Run the simulation with this input file:

```bash
mpirun -np 2 <your_path_to_smilei_installation_folder>/smilei input.py
```

For performance improvement (use of OpenMP parallelization) you may try to export:

```bash
export OMP_NUM_THREADS=2
export OMP_SCHEDULE=dynamic
```

Go back to the directory ``Lecture2``:

```bash
cd ..
```

Check to have all the required packages to use the Python script for plotting:

```bash
pip install numpy scipy matplotlib
```

run the script [plot.py](./plot.py) by typing in the terminal (you may need to use python3 instead of python):

```bash
python plot.py
```

This script prints some plots comparing the results of this code and the results of the 1D Python PIC code in the ``Data`` folder. Now you can visualise the plots in the folder ``Output``.

Now, let's run a 2D simulation:

```bash
cd 2D
```

Run a Smilei simulation with the input file in this folder:

```bash
mpirun -np 2 <your_path_to_smilei_installation_folder>/smilei input.py
```

Then you can plot some figures in the ``PLOTS`` folder:

```bash
cd PLOTS
```

and run all the Python scripts:

```bash
python plot_maps.py 
python plot_scalars.py
python plot_spectra.py 
```

In this folder (accessible from your file explorer), you can see maps of the magnetic field and the electron and ion densities, plots of the evolution of the electromagnetic and kinetic energy, and spectra of electrons at different simulation times.
