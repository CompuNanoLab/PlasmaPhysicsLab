# Particle Orbit Theory

Particle orbit theory studies the dynamics of a single charged particle in assigned electric and magnetic fields. This will be the context to which we apply some of the numerical schemes discussed in the [Introduction to Computational Plasma Physics](./Intro_Comp_Plasma_Phys.md). Charged particle dynamics provides useful hints to understand the physical behaviour of a charged population, and, thus, of a plasma. In the following, we adopt a classical description (not relativistic or quantum) and suppose that the field behaviour is completely given. Let's leave the more accurate self-consistent treatment of electromagnetic fields and particle motion to the numerical studies. The charged particle dynamics comes from solving the following ordinary differential equation:

$$ m \dfrac{d^2\mathbf{x}}{dt^2}=q(\mathbf{E}+\mathbf{v}\times\mathbf{B})$$

In the following, we explore different cases of assigned fields.

## $\mathbf{B}$ uniform and constant

The above differential equation reduces to:

$$ m \dfrac{d^2\mathbf{x}}{dt^2}=m\dfrac{d\mathbf{v}}{dt}=q(\mathbf{v}\times\mathbf{B})$$

We can start to derive useful information on the conserved quantities. Multiplying by $\cdot \mathbf{v}$ and by $\cdot \mathbf{B}$, we get:

$$\dfrac{dv^2}{dt}=0 \quad \mathbf{B}\cdot\dfrac{d\mathbf{v}}{dt} =0$$

Let's define the component of the velocity parallel to the magnetic field as:

$$v_\parallel=\mathbf{v}\cdot\dfrac{\mathbf{B}}{B}$$

$v_\parallel$, $v^2$, and $v_{\perp}^2$ are constant during the particle motion. This suggests to use of cylindrical coordinates. Introducing the **cyclotron frequency** $\omega_c=qB/m$:

$$\dfrac{dv_x}{dt}=\omega_c v_y \quad \dfrac{dv_y}{dt}=-\omega_c v_x \quad v_z=v_\parallel=const$$

which solved finally leads to

$$ x(t)=x_c+\rho_L \sin(\omega t+ \phi_0) \quad y(t)=y_c+\rho_L \cos(\omega t+ \phi_0)$$

where $(x_c,y_c)$ is called guiding center or Larmor centre, and $\rho_L=v_\perp/\omega_c$ is the larmor radius. The solutions describe a circular uniform motion in the $x-y$ plane characterized by the Larmor radius while along the direction of the $\mathbf{B}$ field the motion is uniform. The sign of the charge decides on the direction of the circular motion. The motion is called cyclotron motion. The motion is analogous to a current flowing in a coil and produces a magnetic field opposite to the imposed $\mathbf{B}$ (diamagnetic behaviour).

## $\mathbf{B}$ and $\mathbf{F}$ uniform and constant

In this case, we have an additional force $\mathbf{F}$ of arbitrary orientation that could be the result of a constant electric field applied to the charge. The equation of motion now is:

$$m\dfrac{d\mathbf{v}}{dt}=q(\mathbf{v}\times\mathbf{B})+\mathbf{F}$$

We decompose this equation considering different components:

$$m\dfrac{d\mathbf{v}_\parallel}{dt}=\mathbf{F} _\parallel \quad m\dfrac{d\mathbf{v} _ \perp'}{dt}=q(\mathbf{v} _\perp'\times\mathbf{B}) \quad m\dfrac{d\mathbf{v} _{\perp c}}{dt}=q(\mathbf{v} _{\perp c}\times\mathbf{B})+\mathbf{F}$$

where $\mathbf{v}_\perp=\mathbf{v} _ \perp'+\mathbf{v} _{\perp c}$ so a component follows the circular motion described above while the other is the velocity of the guiding center that now might also move.

We want to find a stationary solution for the guiding center if possible, by imposing that the time derivative of $\mathbf{v} _{\perp c}$ is zero:

$$\mathbf{v} _{\perp c} = \dfrac{\mathbf{F} _\perp \times \mathbf{B}} {q B^2}$$

This equation represents a **drift of the guiding center** orthogonal both to $\mathbf{F} _\perp$ and $\mathbf{B}$. The larmor radius in this case changes based on the position of the particle. If $\mathbf{F}=q\mathbf{E}$ the guiding center drift is called ExB drift and is independent of the particle kind:

$$\mathbf{v} _{\perp c} =\mathbf{v}_E= c \dfrac{\mathbf{E} \times \mathbf{B}} {B^2}$$





















