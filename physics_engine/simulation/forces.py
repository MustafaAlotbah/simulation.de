import numpy as np

from .objects import Particle, Spring


gravitational_constant = 6.67384 * (10 ** -6)  # This constant is NOT realistic!


# Gravitational Forces in Newtonian Mechanics
def apply_gravitational_forces(planet1: Particle, planet2: Particle):
    if planet1.apply_gravitational_forces or planet2.apply_gravitational_forces:
        gravity = - gravitational_constant * planet1.mass * planet2.mass / np.linalg.norm(planet1.position - planet2.position) ** 3
        gravitational_force = gravity * (planet1.position - planet2.position)

        if planet1.apply_gravitational_forces:
            planet1.forces += [gravitational_force]

        if planet2.apply_gravitational_forces:
            planet2.forces += [-gravitational_force]