# Introduction to Computational Plasma Physics

## The plasma physics problem
Most of the field of plasma physics can be encapsulated within Maxwell's equations and the relativistic kinetic equations governing the evolution of the particle distribution function, denoted as $f_a(\mathbf{x},\mathbf{p},t)$, for each distinct plasma species $a$ within a six-dimensional (6D) phase space.

$$\dfrac{\partial f_a}{\partial t}+\dfrac{\mathbf{p}}{\gamma_a m_a} \cdot \dfrac{\partial f_a}{\partial \mathbf{x}}+q_a \left(\mathbf{E}+\dfrac{\mathbf{p}}{\gamma_a m_a}\times\mathbf{B} \right) \cdot \dfrac{\partial f_a}{\partial \mathbf{p}}=C_a$$

$$\nabla \cdot \mathbf{B} = 0$$
     
$$\nabla \times \mathbf{E} =-\dfrac{\partial \mathbf{B}}{\partial t}$$
     
$$\nabla \cdot \mathbf{E}=\dfrac{1}{\epsilon_0} \left(\rho_{ext}+\displaystyle\sum_a q_a \int f_a d\mathbf{p} \right)$$
     
$$\nabla \times \mathbf{B}= \mu_0 \left( \mathbf{J}_{ext}+\displaystyle\sum_a \dfrac{q_a}{m_a} \int \dfrac{\mathbf{p}}{\gamma_a} f_a d\mathbf{p}\right)+\epsilon_0\mu_0 \dfrac{\partial \mathbf{E}}{\partial t}$$ 

The $C_a$ term describes the effect of collisions and $\gamma_{a}=\left(1+p^2/(m_a^2c^2)\right)^{\frac{1}{2}}$. If $C_a$ is set to zero, the system takes the name of the relativistic collisionless Vlasov-Maxwell system. In these equations, the fields dictate how particles move, and the particle motion itself modifies the fields. As a result, this is a system of coupled, nonlinear equations, representing a broad spectrum of physics across various temporal and spatial scales. Over time, many approximations to these equations have been developed, often making them more manageable and allowing one to obtain reasonable results in specific physical situations. Nevertheless, the forefront of theoretical and computational plasma physics continues to strive for a comprehensive kinetic comprehension of plasma dynamics, derived from this complete set of equations. 

In the following sections, we will understand the motivations that lead to the use of computational tools in plasma physics and then quickly survey modern computational techniques employed in this context.

## Why are computational methods central in plasma physics?

We know the fundamental laws (written above) that govern plasma physics, but "we are simply unable to work out their consequences" ([Dawson,1983](https://doi.org/10.1103/revmodphys.55.403)). 

In many fields of study, including plasma physics, experiments can be challenging or impossible, and the concurrent interaction of numerous degrees of freedom makes analytical methods unfeasible. Computer simulations are powerful tools to investigate physical scenarios of this kind. 
One starts constructing a numerical model of the system or theory of interest. Then, one carries out a numerical experiment allowing the system to evolve from some initial conditions following the laws implemented. When required, high-performance computing machines are used.
Computer simulations can give all the desired information on the evolution of the numerical system. These results are then compared with theoretical predictions based on simplified analytic models, with experimental outcomes or with the observation of natural phenomena, or one can use the results to predict the behaviour of unperformed experiments.

Computational plasma physics can provide results of immediate practical interest, for example about the performance of a fusion device, of a plasma-based accelerator, or of an apparatus for generating radiation or processing materials. With simulations, we gain also insight and understanding of fundamental physical aspects, like collective mechanisms of energy and plasma transport across electromagnetic fields, the interaction of the solar wind with planetary magnetospheres, the generation of radiation by energetic plasma, the collapse of a gas cloud to form a star, etc.

So the answer to the starting question is... because some theoretical questions in plasma physics can only be answered with computer simulations which notably provide tools to understand and design experiments, observations, and applications.

## Major computational schemes in plasma physics
In particular, we will focus on the following topics 

- Simulation methods in which the Vlasov-Maxwell (collisionless) equation is solved using the finite-difference-time-domain particle-in-cell method (FDTD PIC).
- Solving fluid equations that are obtained from taking moments of the Vlasov-Maxwell equations and making a closure to approximate the moments not evolved by the fluid equations.
- Directly discretizing the Vlasov-Maxwell equations as a partial differential equation in 6D.
- Monte Carlo methods to

Much of modern computational plasma physics is focused on inventing schemes that preserve at least some of the conservations and other properties of the continuous Vlasov-Maxwell system.
Vlasov codes: This is an area of active research and has applications to the study of turbulence in fusion machines and also exploring fundamental plasma physics in phase-space

## Ordinary Differential Equation Solvers
Particle-in-cell methods are based on pushing macro-particles. These represent the motion of characteristics in phase-space, along which the distribution function is conserved. The macro-particle equations-of-motion are

$$\frac{d\mathbf{x}}{dt} = \mathbf{v}$$ 

$$\frac{d\mathbf{v}}{dt} = \frac{q}{m}(\mathbf{E} + \mathbf{v}\times\mathbf{B})$$

Approach it with finite-difference time domain schemes, ie. derivatives in time are converted to differences that we solve in discrete timesteps.
The most widely used method to solve this system of ODEs is the Boris algorithm.
The Boris algorithm is surprisingly good: it is a second-order, time-centered method that conserves phase-space volume.

other pushers...

### Boris Pusher
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

## Hyperbolic Equations
Hyperbolic equations describe a broad class of physical problems and are essentially characterized by finite propagation speed of disturbances. Examples of hyperbolic equations include Maxwell equations, Euler equations for ideal fluids and ideal MHD equations.

### FDTD
The Yee-cell preserves the underlying geometric structure of Maxwell equations, and ensures that the divergence relations are maintained in the case of vacuum fields. In essence, the electric field is a vector quantity (associated with lines) while the magnetic field is a bi-vector quantity (associated with surfaces). Hence, the most natural representation on a discrete grid utilizes this geometric fact to build a consistent scheme.

### FV
We will focus on finite-volume and discontinuous Galerkin schemes for partial differential equations (PDEs), specifically fluid mechanics (Euler equations) and plasma physics (MHD equations, multi-fluid equations and the Vlasov-Maxwell system). 


## Bibliography

* Tajima, T. (2004). Computational Plasma Physics: With Applications To Fusion And Astrophysics. United Kingdom: Avalon Publishing.
* Kawata, S. (2023). Computational Plasma Science: Physics and Selected Simulation Examples. Germany: Springer Nature Singapore.
* Dawson, J. M. (1983). Particle simulation of plasmas. In Reviews of Modern Physics (Vol. 55, Issue 2, pp. 403–447). American Physical Society (APS). https://doi.org/10.1103/revmodphys.55.403
* Birdsall, C.K., & Langdon, A.B. (1991). Plasma Physics via Computer Simulation (1st ed.). CRC Press. https://doi.org/10.1201/9781315275048
