import numpy as np

from .objects import Particle


earth_gravity = np.array([0, -9.81], dtype='float64')
gravitational_constant = 6.67384 * (10 ** -6)  # This constant is NOT realistic!


# Earth Gravity (On Earth's surface)
def apply_earth_gravity(particle: Particle):
    if particle.apply_earth_gravity:
        particle.velocity += earth_gravity


# Gravitational Forces in Newtonian Mechanics
def apply_gravitational_forces(planet1: Particle, planet2: Particle):
    gravity = - gravitational_constant * planet1.mass * planet2.mass / np.linalg.norm(planet1.position - planet2.position) ** 3
    gravitational_force = gravity * (planet1.position - planet2.position)

    if planet1.apply_gravitational_forces:
        planet1.force += gravitational_force

    if planet2.apply_gravitational_forces:
        planet2.force -= gravitational_force
