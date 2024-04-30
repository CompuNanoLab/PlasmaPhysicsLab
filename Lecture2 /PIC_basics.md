# Particle-In-Cell (PIC) codes

## From the previous lecture..
### Particle-In-Cell (PIC) method

Particle codes are the most successful tool for the simulation of the kinetic dynamics of plasmas. They approximately solve the Vlasov-Maxwell system by replacing (sampling) the distribution function of every plasma species with a collection of computational particles – the macro-particles – each one representing multiple physical particles. Particle-In-Cell (PIC) methods specifically adopt this strategy and make macro-particles evolve self-consistently with electromagnetic fields computed with Maxwell equations on a discrete spatial grid. Since the invention of PIC methods ([Dawson,1983](https://doi.org/10.1103/revmodphys.55.403)), the development of new algorithms and the availability of more powerful computers has allowed continuous progress of PIC simulation from simple, one-dimensional, electrostatic problems to more complex and realistic situations, involving electromagnetic fields in three dimensions and tracking millions of macro-particles.

Their success is manifest in the fact that they are employed in several different contexts:
* laser-plasma acceleration [Fonseca et al. Plasma Physics and Controlled Fusion 50.12 (2008)](https://doi.org/10.1088/0741-3335/50/12/124034)
* nuclear fusion [Yin et al. Physics of Plasmas 16.11 (2009)](https://doi.org/10.1063/1.3250928)
* plasma propulsion ([example](https://gauss-supercomputing.de/))
* astrophysics [Inchingolo et al. 58th Annual Meeting of the APS Division
of Plasma Physics, 8, 61 (2016)](https://loureirogroup.mit.edu/magnetorotational-instability)
* low-temperature plasmas [Tonneau at al. Plasma Sources Science and
Technology 29.11 (2020)](https://doi.org/10.1088/1361-6595/abb3a0)

Usually, PIC codes solve the relativistic Vlasov-Maxwell systems therefore their domain of applicability is relativistic ($T\gg m_ec^2$) and collisionless plasmas in which the time scale of collision events is much larger than the one of plasma oscillations ($\omega_p^{-1} \ll \nu_{coll}^{-1}$). However, collisions, ionization and quantum effects not taken into account by the core PIC algorithm may be included with Monte Carlo strategies.

Let's explore the PIC method in more detail to understand its relation to the Vlasov equation. The starting assumption is that the plasma species distribution function can be approximated in the following way:

$$f_a (\mathbf{x}, \mathbf{p}, t)\approx f_h (\mathbf{x}, \mathbf{p}, t)=\sum_{p=1}^{N_{a}} w_{pa} S(\mathbf{x}- \mathbf{x}_ {pa} (t))\delta(\mathbf{p}-\mathbf{p}_{pa}(t))$$

where $S$ is the **shape function** evenly localized around $\mathbf{x}=0$ and normalized to unity, i.e. its integral over the whole space gives one. Practically, each macro-particle has definite momentum, i.e. is represented by a Dirac delta distribution in $\mathbf{p}$, a finite spatial extension described by the shape function, and a weight $w$ which represents the number of physical particles per macro-particle.

If we put this assumption in the Vlasov equation and calculate its moments we get the equations of motion of the macro-particles. See the [appendix](#appendix) for the detailed calculations.

0<sup>th</sup> moment:

$$ \int \int \left[\dfrac{\partial f_h}{\partial t}+\mathbf{v}\cdot\dfrac{\partial f_h}{\partial \mathbf{x}}+\mathbf{F}_ L\cdot\dfrac{\partial f_ h}{\partial \mathbf{p}}\right] d\mathbf{p}d\mathbf{x}=0 \quad \Longrightarrow \quad \dfrac{dw_{pa}}{dt}=0$$

1<sup>st</sup> moment in $\mathbf{x}$:

$$\int \int \mathbf{x}\left[\dfrac{\partial f_h}{\partial t}+\mathbf{v}\cdot\dfrac{\partial f_h}{\partial \mathbf{x}}+\mathbf{F}_ L\cdot\dfrac{\partial f_ h}{\partial \mathbf{p}}\right] d\mathbf{p}d\mathbf{x}=0 \quad \Longrightarrow \quad\dfrac{ d\mathbf{x} _ {pa}}{dt}=\mathbf{v} _ {pa}$$

2<sup>nd</sup> moment in $\mathbf{p}$:

$$\int \int \mathbf{p}\left[\dfrac{\partial f_h}{\partial t}+\mathbf{v}\cdot\dfrac{\partial f_h}{\partial \mathbf{x}}+\mathbf{F}_ L\cdot\dfrac{\partial f_ h}{\partial \mathbf{p}}\right] d\mathbf{p}d\mathbf{x}=0 \quad \Longrightarrow \quad \dfrac{d \mathbf{p}_ {pa}}{dt}=q_a \left(\mathbf{E}_ {pa}+\mathbf{v}_ {pa}\times\mathbf{B}_ {pa}\right)$$

where:

$$\mathbf{E}_ {pa}(t)=\int \mathbf{E}(\mathbf{x},t)S(\mathbf{x}- \mathbf{x}_ {pa} (t)) d\mathbf{x}$$

$$\mathbf{B}_ {pa}(t)=\int \mathbf{B}(\mathbf{x},t)S(\mathbf{x}- \mathbf{x}_ {pa} (t)) d\mathbf{x}$$

These equations of motion are used inside the PIC algorithm. In the following, we will focus on the simpler scheme of finite-difference-time-domain (FDTD) PIC: we adopt a temporal discretization, a spatial grid in 1D,2D or 3D, and spatial and time derivatives are computed with second-order centred finite differences according to the [**leap-frog scheme**]
(https://en.wikipedia.org/wiki/Leapfrog_integration) that will be discussed later. Starting from initial conditions, every time step $\delta t$, the PIC algorithm follows these steps to evolve the plasma dynamics:

* Evaluation of the current density field on the spatial grid nodes from macro-particle charges, positions and velocities.
* Determination of the electromagnetic field on the grid solving Maxwell equations with the computed current density.
* Interpolation of the fields from the grid nodes to the macro-particles positions.
* Update of macro-particle momenta and positions with the Lorentz force acting on them (particle pushing).

These operations constitute a loop (see **Figure **) that can be reiterated, advancing in time, step by step. Maxwell equations are easily solved by direct integration with finite differences, and it is not necessary to solve divergences equations at each timestep if they are satisfied at the initial time, they remain valid for all times provided that the continuity equation is satisfied for all times. The problem requires boundary conditions for fields and particles. The most used ones are periodic boundary conditions. Particle and field data are saved during the simulation and analysed in the post-processing phase. The following section is devoted to explaining better the last step of the loop devoted to particle pushing.

|![imagine](https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/35ffc328-c75d-4f89-937f-cc4e41cfabf3)|
|:--:| 
|**Figure 2** *Core loop of the particle-in-cell method.*|

#### Particle Pusher

Particle-in-cell methods are based on pushing macro-particles following the macro-particle equations of motion for each species s:

$$\frac{d\mathbf{x_p}}{dt} =  \mathbf{v}_p = \dfrac{\mathbf{p}_p}{\gamma_p m_s}$$ 

$$\frac{d\mathbf{p}_p}{dt} = q_s (\mathbf{E} + \mathbf{v}\times\mathbf{B})$$

Let's approach the problem with finite-difference time-domain schemes, ie. derivatives in time are converted to differences in discrete timesteps. In particular, we use the leap-frog algorithm to discretize in time particle position $\mathbf{x}$ and velocity $\mathbf{v}$ (or momentum $\mathbf{p}$) and then update them at staggered time points so that they "leapfrog" over each other.  These are the discretized equations: 

$$\mathbf{x}^{n+1}=\mathbf{x}^n+ \mathbf{v}^{n+1/2}\Delta t$$  

$$\mathbf{p}^{n+1/2}=\mathbf{p}^{n-1/2}+\mathbf{F}^n\Delta t$$

$$\mathbf{F}^n= q_s (\mathbf{E}^n+\mathbf{v}^n\times\mathbf{B}^n)$$

where the last line reports the Lorentz force $\mathbf{F}^n$ due to the electric and magnetic fields interpolated at the center of the macro-particle of position $\mathbf{x}$ at the time $t^n=n\Delta t$. The velocity, or analogously the momentum, evaluated at integer times appears in the Lorentz force. Since they are known only at half-integer times, an approximation must be used:

$$\mathbf{v}^{n}=\dfrac{\mathbf{u}^n}{\gamma^n} \quad \mathbf{u}^{n}=\dfrac{\mathbf{u}^{n+1/2}+\mathbf{u}^{n-1/2}}{2}$$

Substituting $\mathbf{u}^{n+1/2}$ in the equation for momentum:

$$ \mathbf{u}^{n} =\mathbf{u}^{n-1/2} +\dfrac{q_s}{2m_s} (\mathbf{E}^n+\dfrac{ \mathbf{u}^n}{\gamma^n} \times\mathbf{B}^n) \Delta t $$

This equation for $\mathbf{u}^{n}$ can become explicit introducing some definitions and approximations:

$$\mathbf{b}=\dfrac{q_s \Delta t \mathbf{B}^n}{2m_s\gamma^n}$$

$$\mathbf{u}^-=\mathbf{u}^{n-1/2}+\dfrac{q \Delta t \mathbf{E}^n}{2m}\quad \gamma^n=\sqrt{1+\mathbf{u}^n\cdot\mathbf{u}^n}\approx \sqrt{1+\mathbf{u} ^-\cdot\mathbf{u} ^-}$$

this last approximation is a consequence of the fact that terms proportional to $(\Delta t)^2$ have been ignored consistently with the accuracy in time of the leap-frog scheme. The equation for $\mathbf{u}^{n}$ now is:

$$ \mathbf{u}^{n} =\mathbf{u}^- + \mathbf{u}^{n}\times \mathbf{b} $$

and the final explicit equation can be found with some algebraic passages (hints: apply $\times\mathbf{b}$ and $\cdot\mathbf{b}$ to both sides of this equation to get two useful relations and apply a couple of vector identities):

$$\mathbf{u}^{n}=\dfrac{1}{1+b^2}\left[\mathbf{u} ^ -+\mathbf{u} ^-\times \mathbf{b}+\mathbf{b}(\mathbf{u} ^-\cdot \mathbf{b})\right]$$

where $b$ is the magnitude of $\mathbf{b}$. Now all the ingredients for the advancement of macro-particle position and momentum are given. This scheme is called the **Boris pusher**. This algorithm is a stable, explicit, second-order, and time-centered method that conserves phase-space volume and it is widely used to solve acceleration problems.

## All the other steps...

# Current Deposition

# Field Interpolation

# Maxwell Solver
One of the passages of the PIC loop is the resolution of the Maxwell equations given the current density. It is not necessary to mention the charge density because actually it is not necessary to solve Maxwell divergences equations at each timestep: if they are satisfied at the initial time, they remain valid for all times provided that continuity equation is satisfied for all times. Indeed, using the Maxwell curl equations and a vectorial identity, one can easily verify that:

$$\dfrac{\partial}{\partial t}\left(\nabla \cdot \mathbf{B}\right) =\nabla\cdot\dfrac{\partial\mathbf{B}}{\partial t}=-c\nabla \cdot (\nabla \times \mathbf{E})= 0 $$

$$ \dfrac{\partial}{\partial t}\left(\nabla \cdot \mathbf{E}-4 \pi \rho\right)=\nabla\cdot\dfrac{\partial\mathbf{E}}{\partial t}-4\pi\dfrac{\partial\rho}{\partial t}=\nabla\cdot(c\nabla \times \mathbf{B}-4\pi\mathbf{J})-4\pi\dfrac{\partial\rho}{\partial t}=-4\pi\left(\nabla\cdot\mathbf{J}+\dfrac{\partial\rho}{\partial t}\right)=0
$$

Thus the code will solve only Ampère-Maxwell and Faraday equations. 
In general, there are many numerical methods to solve Maxwell equations, one of the most used method in the context of the PIC codes consists in the use of a second-order finite difference time domain (FDTD)solver on a Yee-lattice, first proposed by K.Y. Yee in \cite{Yee1966}. It gives a solution to the non-approximated differential form of Maxwell equations adopting as numerical approximation spatial and temporal discretization so that derivatives are substituted by finite differences. The solver operates in the temporal domain. 
The spatial and time derivatives are computed with second-order centred finite differences according to what is called the \textit{leap-frog scheme}. More precisely this scheme employs both integer and half-integer points placed at the nodes of a spatial grid and at the nodes of its shifted version respectively. These points can be identified by triplets $(i,j,k)$ and $(i+1/2,j+1/2,+k+1/2)$, with $i$ ranging from $1$ to $N_x$, $j$ from $1$ to $N_y$, and $k$ ranging from $1$ to $N_z$. They are equally spaced by $\Delta x$, $\Delta y$, and  $\Delta z$ respectively along each direction. Time is discretized too, each instant is identified by $t^n$, with $n$ ranging from $1$ to $N_t$ with spacing $\Delta t$ from the next instant. The discretization introduced can be synthesized in this way:

$$ t^n=n\Delta t \quad \quad & t^{n+1/2}=(n+1/2)\Delta t \\\
  x_i=x_{min}+i\Delta x \quad \quad & x_{i+1/2}=x_{min}+(i+1/2)\Delta x \\\
  y_j=y_{min}+j\Delta y \quad \quad & y_{j+1/2}=y_{min}+(j+1/2)\Delta y \\\
  z_k=z_{min}+k\Delta z \quad \quad & z_{k+1/2}=z_{min}+(k+1/2)\Delta z \\\$$

where $min$ values identify the origin of the grid. Therefore, the scheme works on a \textit{staggered grid} consisting of two meshes shifted by half-cell along the three spatial direction; analogous considerations can be made for the one dimensional grid of temporal discretization.

The leap-frog scheme is a second-order method to relate quantities, evaluated on the proposed spatial and temporal grids, with its fist and second derivatives. Calling $a_{ijk}^n$ the starting quantity, $b_{ijk}^n$ its first derivative, and $c_{ijk}^n$ its second derivative, the relationship imposed are the following, according to the type of derivative considered:

$$ \text{Time derivatives} \quad \quad & \text{Spatial derivatives in $x$} \\
  a_{ijk}^{n+1}=a_{ijk}^n+ b_{ijk}^{n+1/2}\Delta t\quad \quad & a^n_{(i+1)jk}=a_{ijk}^n+ b^n_{(i+1/2)jk}\Delta x\\
  b_{ijk}^{n+1/2}=b_{ijk}^{n-1/2}+c_{ijk}^n\Delta t\quad \quad & b^n_{(i+1/2)jk}=b_{(i-1/2)jk}^n+ c^n_{ijk}\Delta x\\
  c_{ijk}^n=C(a_{ijk}^n) \quad \quad & c_{ijk}^n=C(a_{ijk}^n) $$

where $C$ indicates a generic function.
Since the time derivative of the electric field $\mathbf{E}$ is related to the curl of the magnetic field $\mathbf{B}$ and viceversa, it is convenient to approximate the electromagnetic field in such a way that the $\mathbf{B}$ field is shifted in time by $\Delta t/2$ with respect to $\mathbf{E}$ field. Thus the time evolution of $\mathbf{E}$, from step $n$ to step $n + 1$, depends on the values of the $\mathbf{B}$ field at time step $n + 1/2$, while the time evolution of $\mathbf{B}$ from step $n + 1/2$ to step $n + 3/2$ depends on the values of the $\mathbf{E}$ field at time step $n + 1$, and so on. Consequently, the computation of electromagnetic fields goes as:

$$ \mathbf{E}^{n+1}=\mathbf{E}^{n}+\Delta t\left[c(\nabla\times\mathbf{B})^{n+1/2}-4\pi\mathbf{J}^{n+1/2}\right] \quad \mathbf{B}^{n+1/2}=\mathbf{B}^{n-1/2}-c\Delta t(\nabla\times\mathbf{E})^{n} $$

From the spatial point of view, real electromagnetic fields are approximated by the values on the grid shifted by half-cell, while consistently with this choice also $\mathbf{J}$ is approximated at the middle of the cell while $\rho$ is taken on the nodes of the cell. In Figure \ref{fig:Yee} there is a representation of how the field components are arranged on a 3D grid. 
Summarizing, the numerical approximations realized with the leap-frog scheme are:

$$ \mathbf{E}(\mathbf{x},t)\approx\left( {E_x}_ {(i+1/2)jk}^n {E_y}_ {i(j+1/2)k}^n {E_z}_ {ij(k+1/2)}^n\right) \\\
    \mathbf{B}(\mathbf{x},t)\approx\left( {B_x}_ {i(j+1/2)(k+1/2)}^{n+1/2} {B_y}_ {(i+1/2)j(k+1/2)}^{n+1/2} {B_z}_ {(i+1/2)(j+1/2)k}^{n+1/2}\right) \\\
    \mathbf{J}(\mathbf{x},t)\approx\left( {J_x}_ {(i+1/2)jk}^{n+1/2} {J_y}_ {i(j+1/2)k}^{n+1/2} {J_z}_ {ij(k+1/2)}^{n+1/2}\right) \\\
    \rho(\mathbf{x},t)\approx\rho_{ijk}^n $$


![figYee](https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/b2267106-98e3-4aec-b7e0-bcf1e0a02491)
\caption{Representation of the cell of a 3D grid with indications of the positions at which $\mathbf{E}$, $\mathbf{B}$, $\mathbf{J}$ and $\rho$ are computed. Taken from \cite{Derouillat2018}.}

What remains to do is the calculation of the curl of the electric and magnetic field. $(\nabla\times\mathbf{E})^{n}$ must be computed on the same points of the corresponding $\mathbf{B}$ component, while $(\nabla\times\mathbf{B})^{n+1/2}$ must be computed at the points of $\mathbf{E}$. In this calculation spatial derivatives are obtained according to the leap frog scheme:

$$ (\nabla\times\mathbf{E})^{n}=\begin{bmatrix}
    \dfrac{{E_z}_ {i(j+1)(k+1/2)}^n-{E_z}_ {ij(k+1/2)}^n}{\Delta y}-\dfrac{{E_y}_ {i(j+1/2)(k+1)}^n-{E_y}_{i(j+1/2)k}^n}{\Delta z}\\\
    -\dfrac{{E_z}_ {(i+1)j(k+1/2)}^n-{E_z}_ {ij(k+1/2)}^n}{\Delta x}+\dfrac{{E_x}_ {(i+1/2)j(k+1)}^n-{E_x}_ {(i+1/2)jk}^n}{\Delta z}\\\
    \dfrac{{E_y}_ {(i+1)(j+1/2)k}^n-{E_y}_ {i(j+1/2)k}^n}{\Delta x}-\dfrac{{E_x}_ {(i+1/2)(j+1)k}^n-{E_x}_ {(i+1/2)jk}^n}{\Delta y}
    \end{bmatrix}
$$

$$
(\nabla\times\mathbf{B})^{n+1/2}=\begin{bmatrix}
     \dfrac{{B_z}_ {(i+1/2)(j+1/2)k}^{n+1/2}-{B_z}_ {(i+1/2)(j-1/2)k}^{n+1/2}}{\Delta y}-\dfrac{{B_y}_ {(i+1/2)j(k+1/2)}^{n+1/2}-{B_y}_{(i+1/2)j(k-1/2)}^{n+1/2}}{\Delta z}\\\
    -\dfrac{{B_z}_ {(i+1/2)(j+1/2)k}^{n+1/2}-{B_z}_ {(i-1/2)(j+1/2)k}^{n+1/2}}{\Delta x}+\dfrac{{B_x}_ {i(j+1/2)(k+1/2)}^{n+1/2}-{B_x}_{i(j+1/2)(k-1/2)}^{n+1/2}}{\Delta z}\\\
    \dfrac{{B_y}_ {(i+1/2)j(k+1/2)}^{n+1/2}-{B_y}_ {(i-1/2)j(k+1/2)}^{n+1/2}}{\Delta x}-\dfrac{{B_x}_ {i(j+1/2)(k+1/2)}^{n+1/2}-{B_x}_{i(j-1/2)(k+1/2)}^{n+1/2}}{\Delta y}  \end{bmatrix}
$$

In order to guarantee numerical stability, the spatial and time steps must satisfy the \textit{Courant-Friedrichs-Lewy condition}.
The condition is imposed by the fact that if a wave is moving across a discrete spatial grid and we want to compute its amplitude at discrete time steps of equal duration, then this duration must be less than the time for the wave to travel to adjacent grid points. Therefore, the condition is expressed by:

$$
    c \Delta t \sqrt{\dfrac{1}{{\Delta x}^2}+\dfrac{1}{{\Delta y}^2}+\dfrac{1}{{\Delta z}^2}}=CFL<1
$$

where $c$ is the speed of light.

The second-order leap-frog scheme in space and time is easy to implement and efficient, however higher accuracy is gained only significantly increasing the space-time resolution. Alternative schemes exist to get high accuracy with limited computational requirements.

