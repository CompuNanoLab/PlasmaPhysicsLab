# Introduction to Computational Plasma Physics

## The plasma physics problem
Most of the field of plasma physics can be encapsulated within Maxwell's equations and the relativistic kinetic equations governing the evolution of the particle distribution function, denoted as $f_a(\mathbf{x},\mathbf{p},t)$, for each distinct plasma species $a$ within a six-dimensional (6D) phase space.

$$\dfrac{\partial f_a}{\partial t}+\mathbf{v} \cdot \dfrac{\partial f_a}{\partial \mathbf{x}}+ \mathbf{F}_L \cdot \dfrac{\partial f_a}{\partial \mathbf{p}}=C_a$$

$$\nabla \cdot \mathbf{B} = 0$$
     
$$\nabla \times \mathbf{E} =-\dfrac{\partial \mathbf{B}}{\partial t}$$
     
$$\nabla \cdot \mathbf{E}=\dfrac{1}{\epsilon_0} \left(\rho_{ext}+\displaystyle\sum_a q_a \int f_a d\mathbf{p} \right)$$
     
$$\nabla \times \mathbf{B}= \mu_0 \left( \mathbf{J}_{ext}+\displaystyle\sum_a q_a \int \mathbf{v} f_a d\mathbf{p}\right)+\epsilon_0\mu_0 \dfrac{\partial \mathbf{E}}{\partial t}$$ 

The $C_a$ term describes the effect of collisions, $\mathbf{F}_ L = q_a \left(\mathbf{E}+\mathbf{v}\times\mathbf{B} \right)$, and $\mathbf{v}=\mathbf{p}/(\gamma_a m_a)$ with $\gamma_{a}=\left(1+p^2/(m_a^2c^2)\right)^{\frac{1}{2}}$. If $C_a$ is set to zero, we talk about the relativistic collisionless Vlasov-Maxwell system. In these equations, the fields dictate how particles move, and the particle motion itself modifies the fields. As a result, this is a system of coupled, nonlinear equations, representing a broad spectrum of physics across various temporal and spatial scales. Over time, many approximations to these equations have been developed, often making them more manageable and allowing one to obtain reasonable results in specific physical situations. Nevertheless, the forefront of theoretical and computational plasma physics continues to strive for a comprehensive kinetic comprehension of plasma dynamics, derived from this complete set of equations. 

In the following sections, we will understand the motivations that lead to the use of computational tools in plasma physics and then quickly survey modern computational techniques employed in this context.

## Why are computational methods central in plasma physics?

We know the fundamental laws (written above) that govern plasma physics, but "we are simply unable to work out their consequences" ([Dawson,1983](https://doi.org/10.1103/revmodphys.55.403)). 

In many fields of study, including plasma physics, experiments can be challenging or impossible, and the concurrent interaction of numerous degrees of freedom makes analytical methods unfeasible. Computer simulations are powerful tools to investigate physical scenarios of this kind. 
One starts constructing a numerical model of the system or theory of interest. Then, one carries out a numerical experiment allowing the system to evolve from some initial conditions following the laws implemented. When required, high-performance computing machines are used.
Computer simulations can give all the desired information on the evolution of the numerical system. These results are then compared with theoretical predictions based on simplified analytic models, with experimental outcomes or with the observation of natural phenomena, or one can use the results to predict the behaviour of unperformed experiments.

Computational plasma physics can provide results of immediate practical interest, for example about the performance of a fusion device, of a plasma-based accelerator, or of an apparatus for generating radiation or processing materials. With simulations, we gain also insight and understanding of fundamental physical aspects, like collective mechanisms of energy and plasma transport across electromagnetic fields, the interaction of the solar wind with planetary magnetospheres, the generation of radiation by energetic plasma, the collapse of a gas cloud to form a star, etc.

So the answer to the starting question is... because some theoretical questions in plasma physics can only be answered with computer simulations which notably provide tools to understand and design experiments, observations, and applications.

## Major computational schemes in plasma physics

As pointed out, many approximations to the kinetic equations to describe a plasma have been developed to make the problem manageable in specific physical situations. This has led, starting from the full kinetic model to the kinetic collisionless Vlasov-Maxwell system, and the multi-fluid model and magnetohydrodynamics (MHD) where we neglect individual particles in favour of a continuous fluid description. Since plasma exhibits a hierarchical nature with a variety of instabilities and phenomena at various time scales and spatial scales, the adoption of multiple theoretical models suited for different scales seems a reasonable approach. Yet, the problem of understanding plasma collective behaviours remains extremely complex due to the interrelation between phenomena belonging to different layers of this spatial and temporal hierarchical structure.

Computational plasma physics follows the hierarchy of the different plasma models implementing various numerical strategies for each of them. Computer simulations of plasmas comprise two general areas based on kinetic and fluid descriptions connected by hybrid approaches. [Figure 1](\ref{fig1}) is a tentative classification of numerical strategies employed in plasma physics.

|![immagine \label{fig1}](https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/75519b74-0858-42ca-8017-264db1cd09a1)|
|:--:| 
|**Figure 1** *Classification of numerical strategies employed in plasma physics and their properties in relation to plasma physics models.*|

One very accurate procedure consists in solving numerically the plasma kinetic equations. This can be achieved by directly discretizing the Vlasov-Maxwell equations as a partial differential equation in 6D. This is the approach of **Vlasov codes**, an area of active research in computational plasma physics, mainly used to explore fundamental plasma physics in phase space. The Vlasov-Maxwell system can be approximately solved also with **particle codes**. They simply compute the motions of a collection of charged particles, representative of many real particles, interacting with each other and with externally applied fields. The finite-difference-time-domain particle-in-cell method (FDTD PIC) which we will briefly discuss in the following sections belongs to this category. In general, kinetic simulations have been particularly successful in dealing with physical problems in which the particle distributions deviate significantly from a local Maxwellian distribution. However, they ignore collisions and require huge amounts of computational resources limiting their applicability to the whole plasma physics scenario.

**Fluid/MHD codes**, on the other hand, have generally been applied to large-scale problems directly related to the behaviour of experimental devices. Fluid simulations practically solve fluid equations that are obtained from taking moments of the kinetic equations, making a closure to approximate the moments not evolved by the equations themselves. This usually means that some transport coefficients like mobility, diffusion coefficient, and averaged collision frequency must be provided as input parameters. Fluid codes are older and established approaches, conceptually
straightforward and generally exploit the numerical techniques employed in computational fluid dynamics. Besides large scales, fluid codes can be efficiently applied to study dense plasmas, but they fail to capture kinetic effects associated with non-thermal particles.

In between kinetic and fluid codes stand **hybrid codes** that can exploit various strategies. One strategy consists in simultaneously applying fluid and particle treatments to different species of a given plasma, for example considering fluid-like electrons and particle-like ions. In other cases, the hybrid character comes from solving fluid equations with particle methods. In the guiding-center method, fluid-like treatment is adopted in the direction perpendicular to the ambient magnetic field while the dynamics parallel to the magnetic field is particle-like. 

We conclude the survey of computational schemes with a mention of **transport and Monte Carlo methods**. They allow for achieving a numerical description of what happens at long timescales when collisions and diffusion prevail and decide on the transport of energy and particles in the plasma. The transport equations by nature most often take the form of the diffusion equation and can be solved numerically. Monte Carlo methods can be employed in multiple ways and also coupled with the approaches listed above for example to introduce collisional events artificially and to include the creation and annihilation of particles.

## Particle-In-Cell (PIC) method

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

If we put this assumption in the Vlasov equation and calculate its moments we get the equations of motion of the macro-particles.

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
* Interpolation of the fields from the grid nodes to the macro-particles positions;
* Update of macro-particle momenta and positions with the Lorentz force acting on them (particle pushing);

These operations constitute a loop (see [Figure 2](\ref{fig2})) that can be reiterated, advancing in time, step by step. Maxwell equations are easily solved by direct integration with finite differences, and it is not necessary to solve divergences equations at each timestep if they are satisfied at the initial time, they remain valid for all times provided that the continuity equation is satisfied for all times. The problem requires boundary conditions for fields and particles. The most used ones are periodic boundary conditions. Particle and field data are saved during the simulation and analysed in the post-processing phase. The following section is devoted to explaining better the last step of the loop devoted to particle pushing.

### Particle Pusher

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

this last approximation is a consequence of the fact that terms proportional to $(\Delta t)^2$ have been ignored consistently with the accuracy in time of the leap-frog scheme. The equation for $ mathbf{u}^{n} now is:

$$ \mathbf{u}^{n} =\mathbf{u}^- + \mathbf{u}^{n}\times \mathbf{b} $$

and the final explicit equation can be found with some algebraic passages (hints: apply $\times\mathbf{b}$ and $\cdot\mathbf{b}$ to both sides of this equation to get two useful relations and apply a couple of vector identities):

$$\mathbf{u}^{n}=\dfrac{1}{1+b^2}\left[\mathbf{u} ^ -+\mathbf{u} ^-\times \mathbf{b}+\mathbf{b}(\mathbf{u} ^-\cdot \mathbf{b})\right]$$

where $b$ is the magnitude of $\mathbf{b}$. Now all the ingredients for the advancement of macro-particle position and momentum are given. This scheme is called the **Boris pusher**. This algorithm is a stable, explicit, second-order, and time-centered method that conserves phase-space volume and it is widely used to solve acceleration problems.

## Bibliography

* Tajima, T. (2004). Computational Plasma Physics: With Applications To Fusion And Astrophysics. United Kingdom: Avalon Publishing.
* Jardin, S. (2010). Computational Methods in Plasma Physics. United States: CRC Press.
* Dawson, J. M. (1983). Particle simulation of plasmas. In Reviews of Modern Physics (Vol. 55, Issue 2, pp. 403–447). American Physical Society (APS). https://doi.org/10.1103/revmodphys.55.403
* Birdsall, C.K., & Langdon, A.B. (1991). Plasma Physics via Computer Simulation (1st ed.). CRC Press. https://doi.org/10.1201/9781315275048
* [Open source PIC code Smilei](https://smileipic.github.io/Smilei/)
