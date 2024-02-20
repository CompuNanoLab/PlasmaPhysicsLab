# Introduction to Computational Plasma Physics

## The plasma physics problem
Most of the field of plasma physics can be encapsulated within Maxwell's equations and the relativistic kinetic equations governing the evolution of the particle distribution function, denoted as $f_a(\mathbf{x},\mathbf{p},t)$, for each distinct plasma species $a$ within a six-dimensional (6D) phase space.

$$\dfrac{\partial f_a}{\partial t}+\dfrac{\mathbf{p}}{\gamma_a m_a} \cdot \dfrac{\partial f_a}{\partial \mathbf{x}}+q_a \left(\mathbf{E}+\dfrac{\mathbf{p}}{\gamma_a m_a}\times\mathbf{B} \right) \cdot \dfrac{\partial f_a}{\partial \mathbf{p}}=C_a$$

$$\nabla \cdot \mathbf{B} = 0$$
     
$$\nabla \times \mathbf{E} =-\dfrac{\partial \mathbf{B}}{\partial t}$$
     
$$\nabla \cdot \mathbf{E}=\dfrac{1}{\epsilon_0} \left(\rho_{ext}+\displaystyle\sum_a q_a \int f_a d\mathbf{p} \right)$$
     
$$\nabla \times \mathbf{B}= \mu_0 \left( \mathbf{J}_{ext}+\displaystyle\sum_a \dfrac{q_a}{m_a} \int \dfrac{\mathbf{p}}{\gamma_a} f_a d\mathbf{p}\right)+\epsilon_0\mu_0 \dfrac{\partial \mathbf{E}}{\partial t}$$ 

The $C_a$ term describes the effect of collisions and $\gamma_{a}=\left(1+p^2/(m_a^2c^2)\right)^{\frac{1}{2}}$. If $C_a$ is set to zero, we talk about the relativistic collisionless Vlasov-Maxwell system. In these equations, the fields dictate how particles move, and the particle motion itself modifies the fields. As a result, this is a system of coupled, nonlinear equations, representing a broad spectrum of physics across various temporal and spatial scales. Over time, many approximations to these equations have been developed, often making them more manageable and allowing one to obtain reasonable results in specific physical situations. Nevertheless, the forefront of theoretical and computational plasma physics continues to strive for a comprehensive kinetic comprehension of plasma dynamics, derived from this complete set of equations. 

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

Particle codes are the most successful tool for the simulation of the kinetic dynamics of plasmas. 
They approximately solve Vlasov-Maxwell system by replacing (sampling) the distribution function of every species with a
collection of computational particles – the macro-particles – each one
representing multiple physical particles. Particle-In-Cell methods specifically adopt this strategy and make macro-particles evolve self-consistently with electromagnetic fields computed with Maxwell equations on a discrete spatial grid. 
Since the invention of PIC methods (Dawson, 1983), the development of new algorithms and the availability of more powerful computers has allowed continuous progress of PIC simulation from simple, one-dimensional, electrostatic problems to more complex and realistic situations, involving electromagnetic fields in three dimensions and tracking millions of maxro-particles.

Their success is manifest in the fact that they are employed in several different contexts:
* laser-plasma acceleration [Fonseca et al. Plasma Physics and Controlled Fusion 50.12 (2008)](https://doi.org/10.1088/0741-3335/50/12/124034)
* nuclear fusion [Yin et al. Physics of Plasmas 16.11 (2009)](https://doi.org/10.1063/1.3250928)
* plasma propulsion ([example](https://gauss-supercomputing.de/))
* astrophysics [Inchingolo et al. 58th Annual Meeting of the APS Division
of Plasma Physics, 8, 61 (2016)](https://loureirogroup.mit.edu/magnetorotational-instability)
* low-temperature plasmas [Tonneau at al. Plasma Sources Science and
Technology 29.11 (2020)](https://doi.org/10.1088/1361-6595/abb3a0)

Usually, PIC codes solve the relativistic Vlasov-Maxwell systems therefore their domain of applicability are relativistic ($ T>> m_ec^2$) and
collisionless plasmas in which the time scale of collision events is much larger than the one of plasma oscillations ($\omega_p^{-1} \less\less \nu_{coll}^{-1}$). Howeber, collisions, ionization and quantum effects not taken into account by the core PIC algorithm may be included with Monte Carlo strategies.

Let's explore the PIC method in more details to understand it's relation to the Vlasov equation. 

The starting assumption is that the plasma species distribution function can be approximated in the following way

... PIC ansatz...

where ...
Practically, each macro-particle has definite momentum, therefore is represented by a Dirac delta distribution in $\mathbf{p}$,a finite spatial extension described by the shape function $S$ in $\mathbf{x}$, and a weight $w$ which represents the number of physical particles per macro-particle.

If we put this assumption in the Vlasov equation and calculate it's moments we get the equation of motion of the maxro-particles:

equations ....

we use the right-hand side equations inside the PIC algorithm. In the following we will focus on the simpler scheme of finite-difference-time-domain (FDTD) PIC codes having in mind that alternatives exist, for example employing spectral methods, i.e. Fourier transform to solve fields.
So, we simply adopt a temporal and spatial discretization. Every time step $\delta t$ starting from initial conditions, the PIC algorithm follows these steps to evolve the plasma dynamics:

* computes the electric and magnetic fields with Maxwell equations. Maxwell equations are easily solved by direct integration with finite differences on a spatial grid in 1D, 2D or 3D.

* Inter


core loop commented

They are used in 




Much of modern computational plasma physics is focused on inventing schemes that preserve at least some of the conservations and other properties of the continuous Vlasov-Maxwell system.


Hyperbolic equations describe a broad class of physical problems and are essentially characterized by finite propagation speed of disturbances. Examples of hyperbolic equations include Maxwell equations, Euler equations for ideal fluids and ideal MHD equations. The Yee-cell preserves the underlying geometric structure of Maxwell equations, and ensures that the divergence relations are maintained in the case of vacuum fields. In essence, the electric field is a vector quantity (associated with lines) while the magnetic field is a bi-vector quantity (associated with surfaces). Hence, the most natural representation on a discrete grid utilizes this geometric fact to build a consistent scheme.


### Particle Pusher

Particle-in-cell methods are based on pushing macro-particles. These represent the motion of characteristics in phase-space, along which the distribution function is conserved. The macro-particle equations-of-motion are

$$\frac{d\mathbf{x}}{dt} = \mathbf{v}$$ 

$$\frac{d\mathbf{v}}{dt} = \frac{q}{m}(\mathbf{E} + \mathbf{v}\times\mathbf{B})$$

Approach it with finite-difference time-domain schemes, ie. derivatives in time are converted to differences that we solve in discrete timesteps.
The most widely used method to solve this system of ODEs is the Boris algorithm.
The Boris algorithm is surprisingly good: it is a second-order, time-centered method that conserves phase-space volume.

t is an accurate second-
order scheme for the acceleration equation, i.e., the error is proportional to
(At)2 ; (b) It is explicit and very simple (does not require any extra infor-
mation); and (c) It is stable and neutral within the threshold At value.
Leap-frog algorithm is used also for the discretization in time of particle position $\mathbf{x}$ and velocity $\mathbf{v}$ (or momentum $\mathbf{p}$). The integration scheme follows from  considering only temporal discretization, with 

$$\mathbf{x}^{n+1}=\mathbf{x}^n+ \mathbf{v}^{n+1/2}\Delta t$$  

$$\mathbf{p}^{n+1/2}=\mathbf{p}^{n-1/2}+\mathbf{F}^n\Delta t$$

$$\mathbf{F}^n=\dfrac{q}{m}(\mathbf{E}^n+\mathbf{v}^n\times\mathbf{B}^n)$$

where, in the last line, it is written the Lorentz force $\mathbf{F}^n$ due to the electric and magnetic fields interpolated at the centre of the macro-particle of position $\mathbf{x}$ at the time $t^n=n\Delta t$. Here appears also the velocity, or analogously the momentum, evaluated at integer times. Since they are known only at half-integer times, an approximation must be used:

$$\mathbf{v}^{n}=\dfrac{\mathbf{p}^n}{\gamma^n} \quad \mathbf{p}^{n}=\dfrac{\mathbf{p}^{n+1/2}+\mathbf{p}^{n-1/2}}{2}$$

Substituting $\mathbf{p}^{n+1/2}$ in the equation for momentum of, one gets:

$$ \mathbf{p}^{n} =\mathbf{p}^{n-1/2} +\dfrac{q}{2m} (\mathbf{E}^n+\dfrac{ \mathbf{p}^n}{\gamma^n \times\mathbf{B}^n}) \Delta t $$

This equation for $\mathbf{p}^{n}$ can become explicit, introducing some definitions and approximations:

$$\mathbf{b}=\dfrac{q \Delta t \mathbf{B}^n}{2m\gamma^n}$$

$$\Tilde{\mathbf{p}}=\mathbf{p}^{n-1/2}+\dfrac{q \Delta t \mathbf{E}^n}{2m}\quad \gamma^n=\sqrt{1+\mathbf{p}^n\cdot\mathbf{p}^n}\approx \sqrt{1+\Tilde{\mathbf{p}}^n\cdot\Tilde{\mathbf{p}}^n}$$

this last approximation is a consequence of the fact that terms proportional to $(\Delta t)^2$ have been ignored consistently with the accuracy in time of the leap-frog scheme. The explicit equation can be found by multiplying $\times\mathbf{b}$ equation:

$$\mathbf{p}^{n}=\dfrac{1}{1+b^2}\left[\Tilde{\mathbf{p}}+\Tilde{\mathbf{p}}\times \mathbf{b}+\mathbf{b}(\Tilde{\mathbf{p}}\cdot \mathbf{b})\right]$$

where $b$ is the magnitude of $\mathbf{b}$.
Now all the ingredients for the advancement of macro-particle position and momentum are given. This scheme is called \textit{Boris-pusher algorithm}. 

The
Runge-Kutta methods-10 is a popular and powerful method for integrating
nonlinear differential equations.


## Bibliography

* Tajima, T. (2004). Computational Plasma Physics: With Applications To Fusion And Astrophysics. United Kingdom: Avalon Publishing.
* Jardin, S. (2010). Computational Methods in Plasma Physics. United States: CRC Press.
* Dawson, J. M. (1983). Particle simulation of plasmas. In Reviews of Modern Physics (Vol. 55, Issue 2, pp. 403–447). American Physical Society (APS). https://doi.org/10.1103/revmodphys.55.403
* Birdsall, C.K., & Langdon, A.B. (1991). Plasma Physics via Computer Simulation (1st ed.). CRC Press. https://doi.org/10.1201/9781315275048
* [Open source PIC code Smilei](https://smileipic.github.io/Smilei/)
