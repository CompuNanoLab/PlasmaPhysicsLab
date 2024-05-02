import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.constants import micron,c,pi,centi,femto,e,epsilon_0,m_e
import happi 
matplotlib.use('Agg')

#constants
lambda_SI = 0.8*micron # wavelength
omega_SI = 2.0*pi*c / lambda_SI
n_crit = (m_e*epsilon_0*(2*pi*c)**2)/((e*lambda_SI)**2)
my_dpi = 300

#open directory
smilei_dir = '..'
s = happi.Open(smilei_dir)

data = s.Scalar(scalar='Uelm', units=['J/m', 'fs'])
times = (data.getTimes())
Uelm = data.getData()

Ukin_e = s.Scalar(scalar='Ukin_ele', units=['J/m', 'fs']).getData()
Ukin_i = s.Scalar(scalar='Ukin_ion', units=['J/m', 'fs']).getData()
Ukin_tot_part  = np.array(Ukin_i)+np.array(Ukin_e)

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8, 8), dpi=my_dpi, sharex=True)

ax[0][0].plot(times, Uelm)
ax[0][0].set_title('Field Energy')

ax[0][1].plot(times, Ukin_tot_part)
ax[0][1].set_title('Total Kinetic Energy')

ax[1][0].plot(times, Ukin_i)
ax[1][0].set_title('Ion Energy')

ax[1][1].plot(times, Ukin_e)
ax[1][1].set_title('Electron Energy')

for a in ax.reshape(-1):
    a.set_xlabel('Time [fs]')
    a.set_ylabel('Energy [J/m]')
    #a.legend(frameon=False,fancybox=False)

image_file_name='scalars.png' 
plt.tight_layout()
plt.savefig(image_file_name,dpi=my_dpi)
plt.close()


