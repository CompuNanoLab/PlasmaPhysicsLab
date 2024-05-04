import happi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.constants import micron, c, pi, e, epsilon_0, m_e

# Use Agg backend for matplotlib
matplotlib.use('Agg')

# Constants and quantities
wavelength = 0.8e-6  # Wavelength in meters
omega = 2 * pi * c / wavelength  # Laser frequency
nc = (m_e * epsilon_0 * (2 * pi * c)**2) / ((e * wavelength)**2)  # Critical plasma density

unitsim = ['MeV', 's']  # Units for simulation data

# Open Smilei simulation
S = happi.Open('..', reference_angular_frequency_SI=omega)

# Particle binning diagnostic
Diag = S.ParticleBinning(0, units=unitsim)
timesteps = Diag.getAvailableTimesteps()
times = Diag.getTimes()
E = Diag.getAxis('ekin')  # Kinetic energy axis

# Iterate over timesteps
for n, ts in enumerate(timesteps[::10]):
    # Calculate electron energy spectrum
    Wele = Diag.getData(ts)[0] * nc * (S.namelist.Lx * S.namelist.Ly * (c / omega)**2) / (
                m_e * c**2 / e * 1e-6)  # Electron energy spectrum
    fig, ax = plt.subplots(1, 1, figsize=(4, 3.5), dpi=300)
    time = times[n]  # Current time
    ax.plot(E, Wele, label=r'Electrons', color='blue')  # Plot electron spectrum
    ax.set_xlabel(r'$E_p$ [MeV]')  # Set x-axis label
    leg = ax.legend(loc='best', ncol=1, shadow=False, fancybox=False, frameon=False)  # Add legend
    ax.set_yscale(r'log')  # Set y-axis scale to logarithmic
    ax.set_ylabel(r'$dN/dE$ [arb. units]')  # Set y-axis label
    ax.set_title('t = %.2f fs' % time)  # Set title with time
    plt.tight_layout()
    fig.savefig('spectra_%04d.png' % ts)  # Save figure
    plt.close()  # Close the current figure
