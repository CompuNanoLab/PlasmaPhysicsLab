import happi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.constants import micron,c,pi,centi,femto,e,epsilon_0,m_e
matplotlib.use('Agg')

#some quantities
wavelen=0.8e-6
omega=(2*pi*c/wavelen)
nc=(m_e*epsilon_0*(2*pi*c)**2)/((e*wavelen)**2)


unitsim = ['MeV','s']
#open smilei simulation
S=happi.Open('..', reference_angular_frequency_SI=omega)
Diag = S.ParticleBinning(0,units = unitsim)
timesteps= Diag.getAvailableTimesteps()
times=Diag.getTimes()
E=Diag.getAxis('ekin')

for n,ts in enumerate(timesteps[::10]):
    Wele=Diag.getData(ts)[0]*nc*(S.namelist.Lx*S.namelist.Ly*(c/omega)**2)/(m_e*c**2/e*1e-6)
    fig,ax=plt.subplots(1,1,figsize=(4,3.5),dpi=300)
    time=times[n]
    ax.plot(E,Wele, label=r'Electrons',color='blue')
    ax.set_xlabel(r'$E_p$ [MeV]')
    leg = ax.legend(loc='best', ncol=1, shadow=False, fancybox=False,frameon=False)
    ax.set_yscale(r'log')
    ax.set_ylabel(r'$dN/dE$ [arb. units]')
    ax.set_title('t = %.2f fs'% time)
    plt.tight_layout()
    fig.savefig('spectra_%04d.png' % ts)
    plt.close()

