import numpy as np

from .objects import Particle, Spring


earth_gravity = np.array([0, -9.81], dtype='float64')
gravitational_constant = 6.67384 * (10 ** -6)  # This constant is NOT realistic!


# Earth Gravity (On Earth's surface)
def apply_earth_gravity(particle: Particle):
    if particle.apply_earth_gravity:
        particle.force += earth_gravity * particle.mass


# Gravitational Forces in Newtonian Mechanics
def apply_gravitational_forces(planet1: Particle, planet2: Particle):
    if planet1.apply_gravitational_forces or planet2.apply_gravitational_forces:
        gravity = - gravitational_constant * planet1.mass * planet2.mass / np.linalg.norm(planet1.position - planet2.position) ** 3
        gravitational_force = gravity * (planet1.position - planet2.position)

        if planet1.apply_gravitational_forces:
            planet1.force += gravitational_force

        if planet2.apply_gravitational_forces:
            planet2.force -= gravitational_force


def apply_spring_forces(spring: Spring):

    # current length of the spring
    dx = spring.particle1.position - spring.particle2.position
    dl = np.linalg.norm(dx)

    # compute the displacement: ||x1 - x2|| - restLength
    c = dl - spring.rest_length

    # determine the gradient

    if dl > 0.001:
        grad_c = dx / dl

        # spring force
        Fs = -spring.stiffness * c * grad_c

        # damping force
        gradc_wrt_t = np.dot(spring.particle1.velocity - spring.particle2.velocity, grad_c)
        Fd = -spring.damping * gradc_wrt_t * grad_c

        if not spring.particle1.is_fixed:
            spring.particle1.force += Fs + Fd

        if not spring.particle2.is_fixed:
            spring.particle2.force -= Fs + Fd
