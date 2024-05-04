import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.constants import micron, c, pi, femto, e, epsilon_0, m_e
import happi

# Use Agg backend for matplotlib
matplotlib.use('Agg')

# Constants
lambda_SI = 0.8 * micron  # Wavelength in meters
omega_SI = 2.0 * pi * c / lambda_SI  # Laser frequency
n_crit = (m_e * epsilon_0 * (2 * pi * c)**2) / ((e * lambda_SI)**2)  # Critical plasma density
my_dpi = 300  # DPI for figures

# Open directory and retrieve data
smilei_dir = '..' 
s = happi.Open(smilei_dir)  # Open the Smilei output directory

# Extract scalar data
data = s.Scalar(scalar='Uelm', units=['J/m', 'fs'])  # Field energy
times = data.getTimes()  # Time array
Uelm = data.getData()  # Field energy data

# Extract kinetic energy data for electrons and ions
Ukin_e = s.Scalar(scalar='Ukin_ele', units=['J/m', 'fs']).getData()  # Electron kinetic energy
Ukin_i = s.Scalar(scalar='Ukin_ion', units=['J/m', 'fs']).getData()  # Ion kinetic energy
Ukin_tot_part = np.array(Ukin_i) + np.array(Ukin_e)  # Total kinetic energy of particles

# Plotting
fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8, 8), dpi=my_dpi, sharex=True)

# Plot field energy
ax[0][0].plot(times, Uelm)
ax[0][0].set_title('Field Energy')

# Plot total kinetic energy
ax[0][1].plot(times, Ukin_tot_part)
ax[0][1].set_title('Total Kinetic Energy')

# Plot ion kinetic energy
ax[1][0].plot(times, Ukin_i)
ax[1][0].set_title('Ion Energy')

# Plot electron kinetic energy
ax[1][1].plot(times, Ukin_e)
ax[1][1].set_title('Electron Energy')

# Set labels for all subplots
for a in ax.reshape(-1):
    a.set_xlabel('Time [fs]')
    a.set_ylabel('Energy [J/m]')

# Save the figure
image_file_name = 'scalars.png'
plt.tight_layout()
plt.savefig(image_file_name, dpi=my_dpi)
plt.close()
