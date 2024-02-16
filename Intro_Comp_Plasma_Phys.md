# Introduction to Computational Plasma Physics

## The Plasma Physics Problem
Practically all plasma physics is contained in Maxwell equations and the relativistic kinetic equations that describe the evolution of a particle distribution function  $f_a(\mathbf{x},\mathbf{p},t)$ for each plasma species $a$ in a 6D phase-space.

$$\dfrac{\partial f_a}{\partial t}+\dfrac{\mathbf{p}}{\gamma_a m_a} \cdot \dfrac{\partial f_a}{\partial \mathbf{x}}+q_a \left(\mathbf{E}+\dfrac{\mathbf{p}}{\gamma_a m_a}\times\mathbf{B} \right) \cdot \dfrac{\partial f_a}{\partial \mathbf{p}}=C_a$$

$$\nabla \cdot \mathbf{B} = 0$$
     
$$\nabla \times \mathbf{E} =-\dfrac{\partial \mathbf{B}}{\partial t}$$
     
$$\nabla \cdot \mathbf{E}=\dfrac{1}{\epsilon_0} \left(\rho_{ext}+\displaystyle\sum_a q_a \int f_a d\mathbf{p} \right)$$
     
$$\nabla \times \mathbf{B}= \mu_0 \left( \mathbf{J}_{ext}+\displaystyle\sum_a \dfrac{q_a}{m_a} \int \dfrac{\mathbf{p}}{\gamma_a} f_a d\mathbf{p}\right)+\epsilon_0\mu_0 \dfrac{\partial \mathbf{E}}{\partial t}$$ 

The $C_a$ term describes the effect of collisions and $\gamma_{a}=\left(1+p^2/(m_a^2c^2)\right)^{\frac{1}{2}}$. This is a system of coupled (the fields tell the particles how to move, and the particle motion itself modifies the fields), nonlinear equations and describes vast physics that spans an enormous range of temporal and spatial scales. Hence, over the decades many approximations to these equations have been developed that are often more tractable, allowing one to obtain reasonable results in specific physical situations. Still the frontier in  theoretical and computational plasma physics remains the complete kinetic understanding of plasma from these full set of equations.

In this lecture, we will quickly survey modern computational techniques to handle problems in plasma physics. 

- Near first-principles simulations methods in which the Vlasov-Maxwell (collisionless) equation is solved using the finite-difference-time-domain particle-in-cell method (FDTD PIC).
- Solving fluid equations that are obtained from taking moments of the Vlasov-Maxwell equations and making a closure to approximate the moments not evolved by the fluid equations.
- Directly discretizing the Vlasov-Maxwell equations as a partial differential equation in 6D. This is an area of active research and has applications to study of turbulence in fusion machines and also exploring fundamental plasma physics in phase-space.

 Much of modern computational plasma physics is focused on inventing schemes that preserve at least some of the conservations and other properties of the continuous Vlasov-Maxwell system. 

## Why computational methods have become central in plasma physics?

We know the fundamental laws written above that govern plasma physics, but we are simply unable to work out their consequences. Indeed, there are problems for which experiments are difficult or impossible, and the simultaneous interaction of a large number of degrees of freedom makes analytic theoretical treatments impractical. For computer simulation one constructs a numerical model of the system or theory which one wishes to investigate. One then carries out a numerical experiment on a high-speed computer, allowing the system to evolve fromsome initial situation of interest in accordance with the laws used. The computer can give one as much information about the details of the evolution as one desires. One can compare the results of each simulation with theoretical predictions based on simplified analytic models, with experimenta1 observations or with observation of natural phenomena, or one can use the results to predict behavior of unperformed (and often unperformable) experiments.

With computer simulations we can get results of immediate practical interest, for example about the performance of a fusion device, the performance of an accelerator, the performance of an electronic device for generating radiation. Or we gain insight and understanding
of fundamental physical aspects, like Collective mechanisms of energy and plasma transport across a magnetic field, collective mechanisms of transport in a fluid, the nature of hydrodynamic turbulence, the interaction of the solar wind with planetary magnetospheres, the generation of radiation by energetic plasma, the collapse of a gas cloud to form a star, the evolution of a galaxy, and the steps by which a complex chemical reaction takes place.

* Some theoretical questions can only be answered with computer simulations.
* providing tools to understand/design experiments or observations.

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
Substituting $\mathbf{p}^{n+1/2}$ in the equation for momentum of , one gets:

$$\mathbf{p}^{n}=\mathbf{p}^{n-1/2}+\dfrac{q}{2m}(\mathbf{E}^n+\dfrac{\mathbf{p}^n}{\gamma^n \times\mathbf{B}^n) \Delta t$$

This equation for $\mathbf{p}^{n}$ can become explicit, introducing some definitions and approximations:
\begin{equation}
  \mathbf{b}=\dfrac{q \Delta t \mathbf{B}^n}{2m\gamma^n} \quad \Tilde{\mathbf{p}}=\mathbf{p}^{n-1/2}+\dfrac{q \Delta t \mathbf{E}^n}{2m}\quad \gamma^n=\sqrt{1+\mathbf{p}^n\cdot\mathbf{p}^n}\approx \sqrt{1+\Tilde{\mathbf{p}}^n\cdot\Tilde{\mathbf{p}}^n}
\end{equation}
this last approximation is a consequence of the fact that terms proportional to $(\Delta t)^2$ have been ignored consistently with the accuracy in time of the leap-frog scheme. The explicit equation can be found by multiplying $\times\mathbf{b}$ equation \ref{sonfolle3}:
\begin{equation}
  \mathbf{p}^{n}=\dfrac{1}{1+b^2}\left[\Tilde{\mathbf{p}}+\Tilde{\mathbf{p}}\times \mathbf{b}+\mathbf{b}(\Tilde{\mathbf{p}}\cdot \mathbf{b})\right]
\end{equation}
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
