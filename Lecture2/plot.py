import happi
import math
from math import pi,sqrt,sin,tan,log10,exp,cos
import numpy as np
from numpy import random, vectorize, trapz
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import os 
import shutil
import matplotlib.pylab as pl
from matplotlib import rcParams

rcParams.update({'font.size': 14})

matplotlib.use('Agg')

smilei_dir = 'smilei'
python_dir = 'Data'
plot_dir = 'Output'
 
if os.path.exists(plot_dir) is True:
    shutil.rmtree(plot_dir)
    os.mkdir(plot_dir)
else:
    os.mkdir(plot_dir)

s = happi.Open('smilei') 

fwhm = s.namelist.fwhm
dt = s.namelist.dt
dx = s.namelist.dx
shift = s.namelist.shift #int(4*fwhm/dt)
#shift=0
#__________________________________________________________
# ENERGY
   
fig, ax = plt.subplots(1,4,figsize=(18,6), dpi=300, sharex=True) 

myEelm = np.loadtxt(python_dir+'/EM_energy.txt')
myELEkin = np.loadtxt(python_dir+'/kinetic_energy_ELE.txt')
myIONkin = np.loadtxt(python_dir+'/kinetic_energy_ION.txt')

Utot = s.Scalar('Utot') 
Uelm = s.Scalar('Uelm')

Ukin = s.Scalar('Ukin') 
Ukin_ele = s.Scalar('Ukin_ELE')
Ukin_ion = s.Scalar('Ukin_ION')

ax[0].plot(Uelm.getTimes(), Uelm.getData(), lw = 8, label = r'U$_{EM}$ Smilei', color = 'blue')
ax[0].plot(myEelm[:,0],myEelm[:,1], lw = 4, label=r'U$_{EM}$ pythonPIC', color = 'red')

ax[1].plot(Ukin.getTimes(), Ukin.getData(), lw = 8, label = r'U$_{kin}$ Smilei', color = 'blue')
ax[1].plot(myELEkin[:,0],myELEkin[:,1]+myIONkin[:,1], lw = 4, label=r'U$_{kin}$ pythonPIC', color = 'red')

ax[2].plot(Ukin_ele.getTimes(), Ukin_ele.getData(), lw = 8, label = r'U$_{kin} ELE$ Smilei', color = 'blue')
ax[2].plot(myELEkin[:,0],myELEkin[:,1], lw = 4, label=r'U$_{kin} ELE$ pythonPIC', color = 'red')

ax[3].plot(Ukin_ion.getTimes(), Ukin_ion.getData(), lw = 8, label = r'U$_{kin} ION$ Smilei', color = 'blue')
ax[3].plot(myIONkin[:,0],myIONkin[:,1], lw = 4, label=r'U$_{kin} ION$ pythonPIC', color = 'red')
for a in ax.reshape(-1):
    a.set_xlabel('time [code units]')
    a.set_ylabel('energy [code units]') 
    a.legend()
plt.tight_layout(rect=[0,0,1, 1])
plt.savefig(plot_dir+"/ENERGY.png")
plt.close()

#__________________________________________________________
# FIELDS

Ex_s = s.Field(0,'Ex')
Ey_s = s.Field(0,'Ey')
Ez_s = s.Field(0,'Ez') 

Bx_s = s.Field(0,'Bx') 
By_s = s.Field(0,'By') 
Bz_s = s.Field(0,'Bz') 

Jx_s = s.Field(0,'Jx')
Jy_s = s.Field(0,'Jy')
Jz_s = s.Field(0,'Jz') 

Rho_ele = s.Field(0,'Rho_ELE')
Rho_ion = s.Field(0,'Rho_ION')


print(Ex_s.getAvailableTimesteps())
print(shift)
print('*************************************************************')

for ts in Ex_s.getAvailableTimesteps():    
    fig, ax = plt.subplots(4,3,figsize=(40,20), dpi=300, sharex=True)

    myEM = np.loadtxt(python_dir+'/EM_fields_%d.txt' %(ts)) 
    Ex_m = myEM[:,0]
    Ey_m = myEM[:,1]
    Ez_m = myEM[:,2]
    x = np.loadtxt(python_dir+'/grid.txt' %(ts)) 

    ax[0][0].plot(Ex_s.getAxis('x'), Ex_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'blue')
    ax[0][0].plot(x, Ex_m,'-', lw  = 2, label='pythonPIC', color = 'dodgerblue')

    ax[1][0].plot(Ey_s.getAxis('x'), Ey_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'red')
    ax[1][0].plot(x, Ey_m,'-', lw  = 2, label='pythonPIC', color = 'orange')

    ax[2][0].plot(Ez_s.getAxis('x'), Ez_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'green')
    ax[2][0].plot(x, Ez_m,'-', lw  = 2, label='pythonPIC', color = 'limegreen')

    ax[0][0].set_title(r'E$_x$')
    ax[1][0].set_title(r'E$_y$')
    ax[2][0].set_title(r'E$_z$')

    myJ = np.loadtxt(python_dir+'/J_field_%d.txt' %(ts)) 
    Jx_m = myJ[:,0]
    Jy_m = myJ[:,1]
    Jz_m = myJ[:,2]
    x = np.loadtxt(python_dir+'/grid.txt' %(ts)) 

    ax[0][2].plot(Jx_s.getAxis('x'), Jx_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'blue')
    ax[0][2].plot(x, Jx_m,'-', lw  = 2, label='pythonPIC', color = 'dodgerblue')

    ax[1][2].plot(Jy_s.getAxis('x'), Jy_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'red')
    ax[1][2].plot(x, Jy_m,'-', lw  = 2, label='pythonPIC', color = 'orange')

    ax[2][2].plot(Jz_s.getAxis('x'), Jz_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'green')
    ax[2][2].plot(x, Jz_m,'-', lw  = 2, label='pythonPIC', color = 'limegreen')

    ax[0][2].set_title(r'J$_x$')
    ax[1][2].set_title(r'J$_y$')
    ax[2][2].set_title(r'J$_z$')

    Bx_m = myEM[:,3]
    By_m = myEM[:,4]
    Bz_m = myEM[:,5]

    ax[0][1].plot(Bx_s.getAxis('x'), Bx_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'blue')
    ax[0][1].plot(x, Bx_m,'-', lw  = 2, label='pythonPIC', color = 'dodgerblue')

    ax[1][1].plot(By_s.getAxis('x'), By_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'red')
    ax[1][1].plot(x, By_m,'-', lw  = 2, label='pythonPIC', color = 'orange')

    ax[2][1].plot(Bz_s.getAxis('x'), Bz_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'green')
    ax[2][1].plot(x, Bz_m,'-', lw  = 2, label='pythonPIC', color = 'limegreen')
    #ax[1][0].plot(Bz_s.getAxis('x'), Bz_s.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'green')
    #ax[1][0].plot(x, Bz_m,'--', lw  = 2, label='pythonPIC', color = 'limegreen')
    ax[0][1].set_title(r'B$_x$')
    ax[1][1].set_title(r'B$_y$')
    ax[2][1].set_title(r'B$_z$')
    
    myrho = np.loadtxt(python_dir+'/ELE_rho_%d.txt' %(ts)) 
    Rho_ele_m= myrho[:,1]
    myrho = np.loadtxt(python_dir+'/ION_rho_%d.txt' %(ts))
    Rho_ion_m = myrho[:,1]
    x = np.loadtxt(python_dir+'/grid.txt' %(ts)) 

    ax[3][0].plot(Rho_ele.getAxis('x'), Rho_ele.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'blue')
    ax[3][0].plot(x, Rho_ele_m,'-', lw  = 2, label='pythonPIC', color = 'dodgerblue')

    ax[3][1].plot(Rho_ion.getAxis('x'), Rho_ion.getData(ts)[0], lw = 2 , label = 'Smilei', color = 'red')
    ax[3][1].plot(x, Rho_ion_m,'-', lw  = 2, label='pythonPIC', color = 'orange')


    ax[3][0].set_title(r'$\rho_{ELE}$')
    ax[3][1].set_title(r'$\rho_{ION}$')

    for a in ax.reshape(-1):
        a.set_xlabel('x [code units]')
        a.legend() 
        #a.axhline(y=5)
        #a.axhline(y=-5)


    plt.tight_layout(rect=[0,0,1, 1])
    plt.savefig(plot_dir+"/FIELDS_%04d.png" % ts)
    plt.close()
