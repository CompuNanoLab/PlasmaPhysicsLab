import numpy as np
from scipy.constants import c,e,m_e
import matplotlib.pyplot as plt

# function performing the Boris algorithm

def Boris_pusher(q,m,pos0,v0,T,dt,E,B,relativity=True):
    
    # position initialization
    x=pos0[0]
    y=pos0[1]
    z=pos0[2]
    
    # velocity initialization
    vx=v0[0]
    vy=v0[1]
    vz=v0[2]
    
    # definition of timesteps
    timesteps=np.arange(0,int(T/dt)+1)
    
    #definition of times
    times=timesteps*dt

    # variables where to store quantities to save
    xsave=[x]
    ysave=[y]
    zsave=[z]
    vxsave=[vx]
    vysave=[vy]
    vzsave=[vz]
    Exsave=[]
    Eysave=[]
    Ezsave=[]
    Bxsave=[]
    Bysave=[]
    Bzsave=[]
    
    # a useful constant
    const=q*dt/(2*m)
    
    # start iterating over timesteps
    for timestep in timesteps:
        
        # get fields at particle position
        Bp=B(x,y,z,times[timestep])
        Ep=E(x,y,z,times[timestep])
        
        # store fields
        Exsave.append(Ep[0])
        Eysave.append(Ep[1])
        Ezsave.append(Ep[2])
        Bxsave.append(Bp[0])
        Bysave.append(Bp[1])
        Bzsave.append(Bp[2])
        
        v2=vx**2+vy**2+vz**2 # squared velocity at step (n-1/2)
        
        # check if relativistic description is required
        if relativity:
            gamma=1/np.sqrt(1-v2) # Lorentz factor at step (n-1/2)
        else:
            gamma=1
            
        # momentum components u=v*gamma  at step (n-1/2)
        ux=vx*gamma
        uy=vy*gamma
        uz=vz*gamma
        
        # u^- components
        u_x=ux+const*Ep[0]
        u_y=uy+const*Ep[1]
        u_z=uz+const*Ep[2]
        
        # check if relativistic description is required
        if relativity:
            gamma_=np.sqrt(1+u_x**2+u_y**2+u_z**2) # approximated Lorentz factor at step (n)
        else:      
            gamma_=1
        
        # b components
        bx=const*Bp[0]/gamma_
        by=const*Bp[1]/gamma_
        bz=const*Bp[2]/gamma_
        
        b2=bx**2+by**2+bz**2 #b^2
        
        u_timesb=u_x*bx+u_y*by+u_z*bz # u^- x b
        
        # momentum components at step (n)
        ux_n=(u_x+u_y*bz-u_z*by+bx*u_timesb)/(1+b2)
        uy_n=(u_y+u_z*bx-u_x*bz+by*u_timesb)/(1+b2)
        uz_n=(u_z+u_x*by-u_y*bx+bz*u_timesb)/(1+b2)
        
        # velocity components at step (n)
        vx_n=ux_n/gamma_
        vy_n=uy_n/gamma_
        vz_n=uz_n/gamma_
        
        # Lorentz force at step (n)
        Fx=q*(Ep[0]+vy_n*Bp[2]-vz_n*Bp[1])
        Fy=q*(Ep[1]+vz_n*Bp[0]-vx_n*Bp[2])
        Fz=q*(Ep[2]+vx_n*Bp[1]-vy_n*Bp[0])
        
        # momentum components at step (n+1/2)
        ux=ux+Fx*dt/m
        uy=uy+Fy*dt/m
        uz=uz+Fz*dt/m
        
        if relativity:
            gamma=np.sqrt(1+ux**2+uy**2+uz**2) # Lorentz factor at step (n+1/2)
        else:      
            gamma=1
            
        # velocity components at step (n+1/2)
        vx=ux/gamma
        vy=uy/gamma
        vz=uz/gamma
        
        # save the updated velocities
        vxsave.append(vx)
        vysave.append(vy)
        vzsave.append(vz)
        
        # position components at step (n+1)
        x=x+vx*dt
        y=y+vy*dt
        z=z+vz*dt
        
        # save the updated positions
        xsave.append(x)
        ysave.append(y)
        zsave.append(z)
        
    return xsave[:-1],ysave[:-1],zsave[:-1],vxsave[:-1],vysave[:-1],vzsave[:-1],Exsave,Eysave,Ezsave,Bxsave,Bysave,Bzsave,times 

# parameters to define

omega=1 # reference frequency
q=1 # particle charge
m=1 # particle mass

pos0=np.array([0,0,0]) # initial position
v0=np.array([0.1,0,0.1]) # initial velocity

T=200/omega # simulation duration 
dt=(1/omega)/100 # simulation timestep

# magnetic field
def B(x,y,z,t):
    Bx=0
    By=0
    Bz=0.1
    return [Bx,By,Bz]
    
# electric field
def E(x,y,z,t):
    Ex=0
    Ey=0
    Ez=0
    return [Ex,Ey,Ez]
    
# apply the Boris pusher
x,y,z,vx,vy,vz,Ex,Ey,Ez,Bx,By,Bz,t = Boris_pusher(q,m,pos0,v0,T,dt,E,B,relativity=False)

# plot trajectory in space

ax = plt.figure(dpi=300).add_subplot(projection='3d')
ax.plot(x,y,z,color='black',label='trajectory')
ax.scatter(x[0],y[0],z[0],color='red',label='start')
ax.scatter(x[-1],y[-1],z[-1],color='blue',label='end')
ax.set_xlabel('$x$ [-]')
ax.set_ylabel('$y$ [-]')
ax.set_zlabel('$z$ [-]')
plt.legend(loc='upper right',fancybox=False)

