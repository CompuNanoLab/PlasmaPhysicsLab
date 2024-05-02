import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.constants import micron,c,pi,centi,femto,e,epsilon_0,m_e,m_p
import happi
matplotlib.use('Agg')


#constants
lambda_SI = 0.8*micron # wavelength
omega_SI = 2.0*pi*c / lambda_SI
mc2 = 0.510998950e6 
n_crit = (m_e*epsilon_0*(2*pi*c)**2)/((e*lambda_SI)**2)
my_dpi = 300

#open directory
smilei_dir = '..'
s = happi.Open(smilei_dir)

#simulation parameters
um_s = 1/s.namelist.um #conversion to micron
fs_s = 1/omega_SI*1e15 #conversion to fs
Lx = s.namelist.Lx*um_s
Ly = s.namelist.Ly*um_s
unitsim = ["fs", "um", "MeV","J/m"]

# get colormap
ncolors = 256
color_array = plt.get_cmap('seismic')(range(ncolors))
color_array[:,-1] = np.concatenate((np.ones(int(14*ncolors/32)),np.zeros(int(4*ncolors/32)),np.ones(int(14*ncolors/32))))
map_object = LinearSegmentedColormap.from_list(name='diverg',colors=color_array)
matplotlib.colormaps.register(cmap=map_object)
 
data = s.Field(0,'Ex')
steps=data.getTimesteps()

for n,ts in enumerate(steps):
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(10, 5), dpi=my_dpi, sharex=True, sharey=True)
    ne = s.Field(0, "Rho_ele", units = unitsim)
    ni = s.Field(0, "Rho_ion", units = unitsim)
    time=ne.getTimes()[n]
    x=ne.getAxis('x')
    y=ne.getAxis('y')
    rhoele=-ne.getData(ts)[0]
    rhoion=ni.getData(ts)[0]
    bz = s.Field(0,"Bz",units = unitsim)
    Bz=bz.getData(timestep=ts)[0]

    ax[0].set_title('Ele+Bz, t = %.2f fs'% time)
    ax[1].set_title('Ions, t = %.2f fs'% time)

    im=ax[0].imshow(np.flip(rhoele,axis=1).T,cmap='Greys',norm=LogNorm(vmin=1e-3, vmax=rhoele.max()),extent=[x.min(),x.max(),y.min(),y.max()], interpolation='bicubic')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('bottom', size='3%',pad=0)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ [n$_c$]')
    im=ax[0].imshow(np.flip(Bz,axis=1).T,cmap='diverg',extent=[x.min(),x.max(),y.min(),y.max()],vmin=-10, vmax=10, interpolation='bicubic')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('bottom', size='3%',pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')
    
    im=ax[1].imshow(np.flip(rhoion,axis=1).T,cmap='Blues',norm=LogNorm(vmin=1e-3, vmax=rhoion.max()),extent=[x.min(),x.max(),y.min(),y.max()], interpolation='bicubic')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('bottom', size='3%',pad=1.1)
    plt.colorbar(im,cax=cax, orientation='horizontal', label=r'n$_i$ [n$_c$]')
    
    for a in ax.reshape(-1):
        a.set_xlabel(r'x [$\mu$m]')
        a.set_ylabel(r'y [$\mu$m]')

    image_file_name ='maps_%04d.png' % ts
    plt.tight_layout()
    plt.savefig(image_file_name,dpi=my_dpi)
    plt.close()
