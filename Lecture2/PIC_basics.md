# Particle-In-Cell (PIC) codes

## From the previous lecture...
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

These equations of motion are used inside the PIC algorithm. In the following, we will focus on the simpler scheme of finite-difference-time-domain (FDTD) PIC: we adopt a temporal discretization, a spatial grid in 1D, 2D or 3D, and spatial and time derivatives are computed with second-order centred finite differences according to the [**leap-frog scheme**](https://en.wikipedia.org/wiki/Leapfrog_integration) that will be discussed later. Starting from initial conditions, every time step $\delta t$, the PIC algorithm follows these steps to evolve the plasma dynamics:

* Evaluation of the current density field on the spatial grid nodes from macro-particle charges, positions and velocities.
* Determination of the electromagnetic field on the grid solving Maxwell equations with the computed current density.
* Interpolation of the fields from the grid nodes to the macro-particles' positions.
* Update of macro-particle momenta and positions with the Lorentz force acting on them (particle pushing).

These operations constitute a loop (see **Figure 1**) that can be reiterated, advancing in time, step by step. Maxwell equations are easily solved by direct integration with finite differences, and it is not necessary to solve divergences equations at each timestep: if they are satisfied at the initial time, they remain valid for all times provided that the continuity equation is satisfied for all times. The problem requires boundary conditions for fields and particles. The most used ones are periodic boundary conditions. Particle and field data are saved during the simulation and analysed in the post-processing phase. The following section is devoted to explaining better the last step of the loop devoted to particle pushing.

|![imagine](https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/35ffc328-c75d-4f89-937f-cc4e41cfabf3)|
|:--:| 
|**Figure 1** Core loop of the particle-in-cell method.|

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

### Maxwell Solver
One step of the PIC loop involves solving the Maxwell equations given the current density. It's unnecessary to solve Maxwell's divergence equations at each timestep: if they are satisfied at the initial time, they remain valid provided that the continuity equation is satisfied for all times. Indeed, utilizing the Maxwell curl equations and a vector identity, one can readily verify that:

$$\dfrac{\partial}{\partial t}\left(\nabla \cdot \mathbf{B}\right) =\nabla\cdot\dfrac{\partial\mathbf{B}}{\partial t}=-\nabla \cdot (\nabla \times \mathbf{E})= 0 $$

$$ \dfrac{\partial}{\partial t}\left(\nabla \cdot \mathbf{E}-\dfrac{1}{\epsilon_0}  \rho\right)=\nabla\cdot\dfrac{\partial\mathbf{E}}{\partial t}-\dfrac{1}{\epsilon_0}\dfrac{\partial\rho}{\partial t}=$$

$$=\dfrac{1}{\epsilon_0\mu_0} \nabla\cdot(\nabla \times \mathbf{B}-\mu_0\mathbf{J})-\dfrac{1}{\epsilon_0}\dfrac{\partial\rho}{\partial t}=-\dfrac{1}{\epsilon_0}\left(\nabla\cdot\mathbf{J}+\dfrac{\partial\rho}{\partial t}\right)=0 $$

Therefore, the code only needs to solve the Ampère-Maxwell and Faraday equations. There are numerous numerical methods to solve them; we focus on one of the most commonly used methods in PIC codes: a second-order finite difference time domain (FDTD) solver on a Yee lattice. This method provides a solution to the non-approximated differential form of Maxwell's equations by adopting spatial and temporal discretizations, replacing derivatives with finite differences. This solver operates in the temporal domain.

The spatial and temporal derivatives are computed using second-order centred finite differences based on the leap-frog scheme. Specifically, this scheme employs integer points placed at the nodes of a spatial grid and half-integer points placed at the nodes of a shifted spatial grid. These points can be identified by triplets $(i, j, k)$ and $(i+1/2, j+1/2, k+1/2)$, where $i$ ranges from $1$ to $N_x$, $j$ from $1$ to $N_y$, and $k$ from $1$ to $N_z$. They are equally spaced by $\Delta x$, $\Delta y$, and $\Delta z$, respectively, along each direction. Time is also discretized, with each instant identified by $t^n$, where $n$ ranges from $1$ to $N_t$ with a spacing of $\Delta t$ from the next instant. This discretization can be summarized as follows:

$$ t^n=n\Delta t \quad \quad t^{n+1/2}=(n+1/2)\Delta t $$

$$ x_i=x_{min}+i\Delta x \quad \quad x_{i+1/2}=x_{min}+(i+1/2)\Delta x $$

$$ y_j=y_{min}+j\Delta y \quad \quad y_{j+1/2}=y_{min}+(j+1/2)\Delta y $$

$$ z_k=z_{min}+k\Delta z \quad \quad z_{k+1/2}=z_{min}+(k+1/2)\Delta z $$

where $min$ values identify the origin of the grid. Therefore, the scheme works on a **staggered grid** consisting of two meshes shifted by half-cell along the three spatial directions; analogous considerations can be made for the one-dimensional grid of temporal discretization.

The leap-frog scheme relates quantities, evaluated on the proposed spatial and temporal grids, with their first and second derivatives. Calling $a_{ijk}^n$ the starting quantity, $b_{ijk}^n$ its first derivative, and $c_{ijk}^n$ its second derivative, the relationship imposed are the following, according to the type of derivative considered:

$$ \text{Time derivatives in $t$} \quad \quad \text{Spatial derivatives in $x$} $$

$$a_{ijk}^{n+1}=a_{ijk}^n+ b_{ijk}^{n+1/2}\Delta t\quad \quad a^n_{(i+1)jk}=a_{ijk}^n+ b^n_{(i+1/2)jk}\Delta x $$

$$b_{ijk}^{n+1/2}=b_{ijk}^{n-1/2}+c_{ijk}^n\Delta t\quad \quad b^n_{(i+1/2)jk}=b_{(i-1/2)jk}^n+ c^n_{ijk}\Delta x $$

$$ c_{ijk}^n=C(a_{ijk}^n) \quad \quad c_{ijk}^n=C(a_{ijk}^n) $$

where $C$ indicates a generic function.

We can apply this scheme to the curl equations. 

Concerning time, since the time derivative of the electric field $\mathbf{E}$ is related to the curl of the magnetic field $\mathbf{B}$ and vice versa, it is convenient to stagger the electromagnetic field in time. This means that the $\mathbf{B}$ field is shifted in time by $\Delta t/2$ relative to the $\mathbf{E}$ field. Thus, the time evolution of $\mathbf{E}$ from step $n$ to step $n + 1$ depends on the values of the $\mathbf{B}$ field at time step $n + 1/2$, while the time evolution of $\mathbf{B}$ from step $n + 1/2$ to step $n + 3/2$ depends on the values of the $\mathbf{E}$ field at time step $n + 1$, and so on. Consequently, the computation of electromagnetic fields proceeds as follows:

$$ \mathbf{E}^{n+1}=\mathbf{E}^{n}+\dfrac{1}{\epsilon_0\mu_0}\Delta t\left[(\nabla\times\mathbf{B})^{n+1/2}-\mu_0\mathbf{J}^{n+1/2}\right] \quad \mathbf{B}^{n+1/2}=\mathbf{B}^{n-1/2}-\Delta t(\nabla\times\mathbf{E})^{n} $$
Concerning space, electromagnetic fields are approximated by the values on the grid, properly shifted by half a cell, following the relations established by the curl equations. Consequently, $\mathbf{J}$ is also approximated at the middle of the cell, while $\rho$ is taken at the nodes of the cell. **Figure 2** illustrates how the field components can be arranged on a 3D grid.

|<img src="https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/b2267106-98e3-4aec-b7e0-bcf1e0a02491" width="400">|
|:--:| 
|**Figure 2** Representation of a cell in a 3D grid with indications of the positions where $\mathbf{E}$, $\mathbf{B}$, $\mathbf{J}$, and $\rho$ are computed.|

Summarizing, the numerical approximations realized with the leap-frog scheme are:

$$ \mathbf{E}(\mathbf{x},t)\approx\left( {E_x}_ {(i+1/2)jk}^n {E_y}_ {i(j+1/2)k}^n {E_z}_ {ij(k+1/2)}^n\right) $$

$$ \mathbf{B}(\mathbf{x},t)\approx\left( {B_x}_ {i(j+1/2)(k+1/2)}^{n+1/2} {B_y}_ {(i+1/2)j(k+1/2)}^{n+1/2} {B_z}_ {(i+1/2)(j+1/2)k}^{n+1/2}\right) $$
    
$$ \mathbf{J}(\mathbf{x},t)\approx\left( {J_x}_ {(i+1/2)jk}^{n+1/2} {J_y}_ {i(j+1/2)k}^{n+1/2} {J_z}_ {ij(k+1/2)}^{n+1/2}\right) $$
    
$$ \rho(\mathbf{x},t)\approx\rho_{ijk}^n $$

Lastly, we have to consider the calculation of the curl of the electric and magnetic fields. $(\nabla\times\mathbf{E})^{n}$ must be computed on the same points of the corresponding $\mathbf{B}$ component, while $(\nabla\times\mathbf{B})^{n+1/2}$ must be computed at the points of $\mathbf{E}$. In this calculation spatial derivatives are obtained according to the leapfrog scheme:

$$ (\nabla\times\mathbf{E})^{n}=\begin{bmatrix}
    \dfrac{{E_z}_ {i(j+1)(k+1/2)}^n-{E_z}_ {ij(k+1/2)}^n}{\Delta y}-\dfrac{{E_y}_ {i(j+1/2)(k+1)}^n-{E_y}_ {i(j+1/2)k}^n}{\Delta z}\\\
    -\dfrac{{E_z}_ {(i+1)j(k+1/2)}^n-{E_z}_ {ij(k+1/2)}^n}{\Delta x}+\dfrac{{E_x}_ {(i+1/2)j(k+1)}^n-{E_x}_ {(i+1/2)jk}^n}{\Delta z}\\\
    \dfrac{{E_y}_ {(i+1)(j+1/2)k}^n-{E_y}_ {i(j+1/2)k}^n}{\Delta x}-\dfrac{{E_x}_ {(i+1/2)(j+1)k}^n-{E_x}_ {(i+1/2)jk}^n}{\Delta y}
    \end{bmatrix} $$

$$ (\nabla\times\mathbf{B})^{n+1/2}=\begin{bmatrix}
     \dfrac{{B_z}_ {(i+1/2)(j+1/2)k}^{n+1/2}-{B_z}_ {(i+1/2)(j-1/2)k}^{n+1/2}}{\Delta y}-\dfrac{{B_y}_ {(i+1/2)j(k+1/2)}^{n+1/2}-{B_y}_ {(i+1/2)j(k-1/2)}^{n+1/2}}{\Delta z}\\\
    -\dfrac{{B_z}_ {(i+1/2)(j+1/2)k}^{n+1/2}-{B_z}_ {(i-1/2)(j+1/2)k}^{n+1/2}}{\Delta x}+\dfrac{{B_x}_ {i(j+1/2)(k+1/2)}^{n+1/2}-{B_x}_ {i(j+1/2)(k-1/2)}^{n+1/2}}{\Delta z}\\\
    \dfrac{{B_y}_ {(i+1/2)j(k+1/2)}^{n+1/2}-{B_y}_ {(i-1/2)j(k+1/2)}^{n+1/2}}{\Delta x}-\dfrac{{B_x}_ {i(j+1/2)(k+1/2)}^{n+1/2}-{B_x}_ {i(j-1/2)(k+1/2)}^{n+1/2}}{\Delta y}  \end{bmatrix} $$

To ensure numerical stability, the spatial and temporal steps must satisfy the **Courant-Friedrichs-Lewy condition**.
This condition is imposed by the requirement that if a wave is traversing a discrete spatial grid and we wish to compute its amplitude at discrete time steps of equal duration, this duration must be shorter than the time it takes for the wave to travel to adjacent grid points. Thus, the condition is expressed as:

$$ c \Delta t \sqrt{\dfrac{1}{{\Delta x}^2}+\dfrac{1}{{\Delta y}^2}+\dfrac{1}{{\Delta z}^2}}=CFL<1 $$

in 1D is:

$$ c \Delta t < \Delta x $$

where $c$ is the speed of light.

### Field Gathering
In the PIC algorithm, the electromagnetic field is computed on the grid, whereas all the quantities related to the macro-particles are evaluated at the particle positions. To compute the Lorentz force acting on each macro-particle at time $t^n$, electric and magnetic fields must be interpolated at the particle position and evaluated at the correct time. The electric field is already known at time $t^n$ so we need only to interpolate it at macro-particle positions $\mathbf{x}_p^n$:
 
$$ \mathbf{E}^n(\mathbf{x}_ p)=\int \mathbf{E}^n S(\mathbf{x}-\mathbf{x}_ {p}^n(t)) d\mathbf{x} $$
 
The magnetic field is known only at half-integer times so it must be first evaluated at the correct time with an approximation and then it can be interpolated as above:

$$ \mathbf{B}^n=\dfrac{\mathbf{B}^{n+1/2}+\mathbf{B}^{n-1/2}}{2} $$

$$ \mathbf{B}^n(\mathbf{x}_ p)=\int \mathbf{B}^n S(\mathbf{x}-\mathbf{x}_{p}^n(t)) d\mathbf{x} $$

The interpolation is numerically performed by combining the field values at the nearest grid points to the particle. This combination depends on the specific shape function $S$. A common choice for shape functions is B-splines. A spline function $B_n(x)$ of order $n$ is a piecewise polynomial function of degree $(n-1)$ in the variable $x$. The points where the pieces meet are referred to as knots, and in the PIC implementation, the spacing between knots matches the grid spacing. In **Figure 3**, B-splines up to third order are depicted in one dimension. These functions are normalized to unity and have units of inverse volume. The shape function of second-order interpolation is typically used to represent macro-particles, as the first-order shape function has a too-sharp profile, while the third-order one has tails that are too long.

|<img src="https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/4b6ca82a-c96f-4613-865f-c7dc2fe11a70" width="400">|
|:--:| 
|**Figure 3** B-splines $B_n(x)$ with $n=0,1,2,3$.|

During the procedure of field gathering, using second-order splines, the contributions to the field felt by the particle come from the nearest grid point and the next adjacent points. These contributions are weighted proportionally to the volume of the shape function contained in the cell around each grid point.

For example, if the particle sits exactly at a grid point in 1D, the contributions to the field felt by the particle come for 75% from the nearest point and 12.5% from each of the other two adjacent points. 

### Current Deposition
The final step to discuss involves reconstructing the charge and current densities on the grid from the velocities and positions of the particles. The current density is necessary to advance the Maxwell solver.
A straightforward method for projecting the current onto the grid involves summing the weighted contribution to the current of each macro-particle at a given grid point. Essentially, a portion of the particle's current is allocated to the nearest grid points, as illustrated in **Figure 4**.

|<img src="https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/8d63ba45-7e4d-48c4-9864-8608d4432b4a" width="400">|
|:--:| 
|**Figure 4** 1D scheme of current and charge deposition.|

Although this procedure is used to calculate charge density and is consistent with the single-particle energy balance, it does not guarantee the continuity equation, i.e. charge conservation, because of unavoidable discrepancies between the weighting on the grid of the charge and current density fields. A charge-conserving algorithm for the deposition of the current density vector from the particles to the grid can be obtained by forcing the continuity equation. Local current density is computed from the discrete version of the continuity equation for the single macro-particle, and then the contributions of all macro-particles are summed on each grid point. This scheme is known as **Esirkepov method**.

### Boundary and Initial Conditions
To conclude the overview of the PIC method, it's important to mention the conditions required for initializing particles and fields to start the simulation and to manage the boundaries of the simulation grid. Fields can be initialized with static contributions or through laser injection. Typical boundary conditions for fields include injecting/absorbing conditions like the Silver-Muller ones, periodic conditions, or reflective conditions. On the other hand, particle initialization involves setting initial positions, momenta, and weights for all macro-particles. At boundary conditions, particles can be removed, reflected, absorbed, thermalized, or reinjected by imposing periodicity.

## <a name="appendix"></a>Appendix

Let's write explicitly the moment of zero order in the $x-y-z$ components:

$$ \int \int \left[\dfrac{\partial f_h}{\partial t}+v_x\dfrac{\partial f_h}{\partial x}+v_y\dfrac{\partial f_h}{\partial y}+v_z\dfrac{\partial f_h}{\partial z}+F_ {Lx}\dfrac{\partial f_ h}{\partial px}+F_ {Ly}\dfrac{\partial f_ h}{\partial py}+F_ {Lz}\dfrac{\partial f_ h}{\partial pz}\right] dp_xdp_ydp_zdxdydz=0$$

$\int\frac{\partial f_h}{\partial x}dx$, $\int\frac{\partial f_h}{\partial p_x}dp_x$,  and analogous terms in $y$ and $z$ are zero because we require a compact support for the shape function $S$ like the Dirac delta (i.e. it is zero outside a small range). The conservation of the macro-particle weight $w_{pa}$ follows from the first term where the time derivative can be taken out of the integrations.

In the moment of first order, due to the same reasons, the only non-zero elements are the element

$$\dfrac{\partial}{\partial t} \int \int \mathbf{x} f_h d\mathbf{p}d\mathbf{x}= \sum w_{pa} \dfrac{\partial}{\partial t} \int \mathbf{x} S(\mathbf{x}- \mathbf{x}_ {pa} (t)) d\mathbf{x}=$$

$$=\sum w_{pa} \left( \dfrac{\partial}{\partial t} \int (\mathbf{x}- \mathbf{x}_ {pa} (t)) S(\mathbf{x}- \mathbf{x}_ {pa} (t)) d\mathbf{x} + \dfrac{\partial}{\partial t} \int \mathbf{x}_ {pa} (t) S(\mathbf{x}- \mathbf{x}_ {pa} (t)) d\mathbf{x} \right)= \sum w_{pa} \dfrac{\partial\mathbf{x}_ {pa} (t)}{\partial t} $$

where the symmetry of $S$ has been used to delete one term, and elements like the following for every coordinate $x$, $y$ and $z$:

$$\int \int x v_x \dfrac{\partial f_h}{\partial x} d\mathbf{p}d\mathbf{x}=-\int \int  v_x  f_h d\mathbf{p}d\mathbf{x} = - \sum w_{pa} v_{pax}$$

using integration by part in $x$. Considering all components, these terms lead to the first equation of motion for the macroparticle. The second-order moment has the term

$$\dfrac{\partial}{\partial t} \int \int \mathbf{p} f_h d\mathbf{p}d\mathbf{x}= \sum w_{pa} \dfrac{\partial\mathbf{p}_ {pa} (t)}{\partial t} $$

using the property of the Dirac delta inside $f_h$, and other non-zero terms like:

$$\int \int p_x F_{Lx} \dfrac{\partial f_h}{\partial p_x} d\mathbf{p}d\mathbf{x}=-\int \int  F_{Lx}  f_h d\mathbf{p}d\mathbf{x} = - \sum w_{pa} q_a \int\left(E_x(\mathbf{x},t)+(\mathbf{v}_ {pa}\times\mathbf{B}(\mathbf{x},t))_ x\right) S(\mathbf{x}- \mathbf{x}_ {pa} (t)) d\mathbf{x}$$

exploiting integration by parts. The second equation of motion is easily retrieved considering the same relations for all components.

## Bibliography

* Dawson, J. M. (1983). Particle simulation of plasmas. In Reviews of Modern Physics (Vol. 55, Issue 2, pp. 403–447). American Physical Society (APS). https://doi.org/10.1103/revmodphys.55.403
* Birdsall, C.K., & Langdon, A.B. (1991). Plasma Physics via Computer Simulation (1st ed.). CRC Press. https://doi.org/10.1201/9781315275048
* [Open source PIC code Smilei](https://smileipic.github.io/Smilei/)
* Esirkepov, T. Zh. (2001). Exact charge conservation scheme for Particle-in-Cell simulation with an arbitrary form-factor. In Computer Physics Communications (Vol. 135, Issue 2, pp. 144–153). Elsevier BV. [https://doi.org/10.1016/s0010-4655(00)00228-9](https://doi.org/10.1016/s0010-4655(00)00228-9)
