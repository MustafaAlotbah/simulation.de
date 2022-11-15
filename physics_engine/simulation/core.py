"""
    The core module for the simulation for simulation.de visualizations
    Mustafa Alotbah 2022
"""

import numpy as np

from .objects import Particle
from .forces import (
    apply_gravitational_forces,
    apply_earth_gravity,
    apply_spring_forces
)


# time intergration method
def apply_symplectic_euler(particle: Particle, dt: float):
    force = np.sum(particle.forces, axis=0)
    particle.velocity += dt * force / particle.mass
    particle.position += dt * particle.velocity


def apply_exlicit_euler(particle: Particle, dt: float):
    force = np.sum(particle.forces, axis=0)
    particle.position += dt * particle.velocity
    particle.velocity += dt * force / particle.mass


# Simulation class
class Simulation:
    def __init__(self, dt=0.1, callback=None, stats_callback=None):
        self.particles = []
        self.springs = []
        self.dt = dt
        self.time = 0.0
        self.callback = callback
        self.stats_callback = stats_callback

    def get_stats(self):
        return str(
            "Stats:\n"
            f"$t={self.time:.2f}$\n"
        ) + (self.stats_callback(self) if self.stats_callback else "")

    def simulation_step(self):
        # calculate forces

        # 1. reset forces
        for particle in self.particles:
            # particle.force = np.array([0.0, 0.0])
            particle.forces = []

        for spring in self.springs:
            apply_spring_forces(spring)

        if self.callback:
            self.callback(self)

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
        for spring in self.springs:
            artists += spring.get_artists()
        for particle in self.particles:
            artists += particle.get_artists()
        return artists

    def get_legend_handles(self):
        artists = []
        artists += [particle.get_legend_handle() for particle in self.particles if particle.get_legend_handle()]
        return artists

    def update_artists(self):
        for particle in self.particles:
            particle.update_artist()
        for spring in self.springs:
            spring.update_artist()

    def get_force_quivers(self):
        offsets, us, vs = [], [], []
        for particle in self.particles:
            p_offsets, p_us, p_vs = particle.get_forces_quivers()
            offsets += p_offsets
            us += p_us
            vs += p_vs
        return offsets, us, vs

