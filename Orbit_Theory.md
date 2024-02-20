# Particle Orbit Theory

Particle orbit theory studies the dynamics of a single charged particle in assigned electric and magnetic fields. This will be the context to which we apply some of the numerical schemes discussed in the [Introduction to Computational Plasma Physics](./Intro_Comp_Plasma_Phys.md). Charged particle dynamics provides useful hints to understand the physical behaviour of a charged population, and, thus, of a plasma. In the following, we adopt a classical description (not relativistic or quantum) and suppose that the field behaviour is completely given. Let's leave the more accurate self-consistent treatment of electromagnetic fields and particle motion to the numerical studies. The charged particle dynamics comes from solving the following ordinary differential equation:

$$ m \dfrac{d^2\mathbf{x}}{dt^2}=q(\mathbf{E}+\mathbf{v}\times\mathbf{B})$$

In the following, we explore different cases of assigned fields.

## $\mtahbf{B}$ uniform and constant

The above differential equation reduces to:

$$ m \dfrac{d^2\mathbf{x}}{dt^2}=m\dfrac{d\mathbf{v}}{dt}=q\mathbf{v}\times\mathbf{B}$$

Instead of solving it (we will do it numerically), we can try to derive useful info on the conserved quantities. Multiplying by $\cdot \mathbf{v}$ and by $\cdot \mathbf{B}$, we get:

$$\dfrac{dv^2}{dt}=0 \quad \mathbf{B}\cdot\dfrac{d\mathbf{v}}{dt} =0$$

Let's define the component of the velocity parallel to the magnetic field as:

$$v_{//}=\mathbf{v}\cdot\dfrac{\mathbf{B}}{B}$$

$v_{//}$, $v^2$, and $v_{\perp}^2$ are constant during the particle motion,







