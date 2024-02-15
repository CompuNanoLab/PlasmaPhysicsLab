# Introduction to Computational Plasma Physics
## 
Vast majority of plasma physics is contained in the Vlasov-Maxwell equations that describes the evolution of a particle distribution $f_s(t,\mathbf{x},\mathbf{v})$ function in 6D phase-space.

$$\dfrac{\partial f_s}{\partial t}+ \nabla_\mathbf{x} \cdot (\mathbf{v}f_s)+ \nabla_\mathbf{v} \cdot \left(\dfrac{q_s}{m_s}(\mathbf{E}+\mathbf{v}\times\mathbf{B}\right) f_s)=\left( \dfrac{\partial f_s}{\partial t} \right)_c$$

where the right-hand side describes the effect of collisions. Of course, the electromagnetic fields are determined from Maxwell equations:

$$\begin{split}\epsilon_0\mu_0 \frac{\partial \mathbf{E}}{\partial t}- \nabla_\mathbf{x} \times \mathbf{B} &= -\mu_0  \sum_s q_s \int_{-\infty}^{\infty} v f_s d\mathbf{v}^3 \\ \frac{\partial \mathbf{B}}{\partial t}+ \nabla_\mathbf{x} \times \mathbf{E} &= 0 \\ \nabla_\mathbf{x}\cdot\mathbf{E} &=\frac{1}{\epsilon_0}\sum_s q_s \int_{-\infty}^{\infty} f_s d\mathbf{v}^3 \\ \nabla_\mathbf{x}\cdot\mathbf{B} &= 0.\end{split}$$

Suitable modifications are required to account for relativistic effects. the problem is highly nonlinear: the fields tell the particles how to move, and the particle motion itself modifies the fields. The Vlasov-Maxwell system is a formidable system of coupled, nonlinear equations and describes vast physics that spans an enormous range of temporal and spatial scales. Hence, over the decades many approximations to the Vlasov-Maxwell equations have been developed that are often more tractable, allowing one to obtain reasonable results in specific physical situations. Of course, the frontier in computational and theoretical plasma physics remains the complete kinetic understanding of plasma from the full Vlasov-Maxwell equations.

In this lecture, we will quickly look at the key modern computational techniques to handle various problems in plasma physics. 

- Near first-principles simulations methods in which the Vlasov-Maxwell equation is solved using the finite-difference-time-domain particle-in-cell method (FDTD PIC).
- Solving fluid equations that are obtained from taking moments of the Vlasov-Maxwell equations and making a closure to approximate the moments not evolved by the fluid equations. Here we will study the modern approach based on the theory of hyperbolic PDEs and Riemann solvers. We will also look at the special requirements of solvers that are required to study fusion problems.
- Directly discretizing the Vlasov-Maxwell equations as a PDE in 6D. This is an emerging area of active research and has applications to study of turbulence in fusion machines and also exploring fundamental plasma physics in phase-space.

 Much of modern computational plasma physics is focused on inventing schemes that preserve at least some of the conservations and other properties of the continuous Vlasov-Maxwell system. 

** Why computational methods have become central in plasma physics? **

We know the fundamental laws written above that govern plasma physics, but we are simply unable to work out their consequences. Indeed, there are problems for which experiments are difficult or impossible, and the simultaneous interaction of a large number of degrees of freedom makes analytic theoretical treatments impractical. For computer simulation one constructs a numerical model of the system or theory which one wishes to investigate. One then carries out a numerical experiment on a high-speed computer, allowing the system to evolve fromsome initial situation of interest in accordance with the laws used. The computer can give one as much information about the details of the evolution as one desires. One can compare the results of each simulation with theoretical predictions based on simplified analytic models, with experimenta1 observations or with observation of natural phenomena, or one can use the results to predict behavior of unperformed (and often unperformable) experiments.

With computer simulations we can get results of immediate practical interest, for example about the performance of a fusion device, the performance of an accelerator, the performance of an electronic device for generating radiation. Or we gain insight and understanding
of fundamental physical aspects, like Collective mechanisms of energy and plasma transport across a magnetic field, collective mechanisms of transport in a fluid, the nature of hydrodynamic turbulence, the interaction of the solar wind with planetary magnetospheres, the generation of radiation by energetic plasma, the collapse of a gas cloud to form a star, the evolution of a galaxy, and the steps by which a complex chemical reaction takes place.



* Some theoretical questions can only be answered with computer simulations.
* providing tools to understand/design experiments or observations.
* 

## Bibliography

* Tajima, T. (2004). Computational Plasma Physics: With Applications To Fusion And Astrophysics. United Kingdom: Avalon Publishing.
* Kawata, S. (2023). Computational Plasma Science: Physics and Selected Simulation Examples. Germany: Springer Nature Singapore.
* Dawson, J. M. (1983). Particle simulation of plasmas. In Reviews of Modern Physics (Vol. 55, Issue 2, pp. 403–447). American Physical Society (APS). https://doi.org/10.1103/revmodphys.55.403
* Birdsall, C.K., & Langdon, A.B. (1991). Plasma Physics via Computer Simulation (1st ed.). CRC Press. https://doi.org/10.1201/9781315275048
* [Link to an online course by Ammar Hakim](https://cmpp.readthedocs.io/en/latest/)
