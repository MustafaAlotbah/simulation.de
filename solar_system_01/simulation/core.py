"""
    The core module for the simulation for simulation.de visualizations
    Mustafa Alotbah 2022
"""

import numpy as np

from .objects import Particle
from .forces import apply_gravitational_forces, apply_earth_gravity


# time intergration method
def apply_symplectic_euler(particle: Particle, dt: float):
    particle.velocity += dt * particle.force / particle.mass
    particle.position += dt * particle.velocity


# Simulation class
class Simulation:
    def __init__(self, dt=0.1):
        self.particles = []
        self.dt = dt
        self.time = 0.0

    def simulation_step(self):
        # calculate forces

        # 1. reset forces
        for particle in self.particles:
            particle.force = np.array([0.0, 0.0])

        for i in range(len(self.particles)):
            # add earth gravity to all particles
            apply_earth_gravity(self.particles[i])

            # add gravitational forces to all unique pairs
            for j in range(len(self.particles)):
                if i < j:
                    apply_gravitational_forces(self.particles[i], self.particles[j])

            # apply symplectic Euler's equations for time integration
            apply_symplectic_euler(self.particles[i], self.dt)

        self.time += self.dt

    #
    def get_artists(self):
        artists = []
        for particle in self.particles:
            artists += particle.get_artists()
        return artists

    def get_legend_handles(self):
        artists = []
        artists += [particle.get_legend_handle() for particle in self.particles]
        return artists

    def update_artists(self):
        for particle in self.particles:
            particle.update_artist()

