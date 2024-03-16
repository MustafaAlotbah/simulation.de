"""
    The core module for the Simulation for simulation.de visualizations.
    This module provides the necessary classes and functions to set up and run physics simulations,
    allowing for easy visualization and analysis of different physical systems.
    Author: Mustafa Alotbah, 2024
"""
import time

import numpy as np
from typing import *
from matplotlib.artist import Artist

from .objects import Particle, Spring
from .forces import apply_gravitational_forces

# Various integrator methods available for the simulation
INTEGRATOR_EXPLICIT_EULER = "explicit-euler"
INTEGRATOR_SYMPLECTIC_EULER = "symplectic-euler"
INTEGRATOR_LEAPFROG = "leapfrog"
INTEGRATOR_RK2 = "runge-kutta-2"
INTEGRATOR_RK4 = "runge-kutta-4"


# The model can also directly use accelerations instead of forces
def compute_acceleration(particle: Particle) -> float:
    """
    Computes the acceleration of a given particle.

    This function takes a Particle object and calculates its acceleration based on the
    net force acting on it and its mass.

    Parameters:
    particle (Particle):    The particle for which acceleration is calculated.

    Returns:
    float:                  The acceleration vector of the particle.
    """
    force = np.sum(particle.forces, axis=0)
    acceleration = force / particle.mass
    return acceleration


class SimulationState:
    """
    Represents the state of the simulation at any given time.

    This class holds all the particles and springs involved in the simulation,
    along with methods to apply various forces and update the state based on
    different integration methods.

    Attributes:
    particles (List[Particle]):     A list of particles in the simulation.
    springs (List[Spring]):         A list of springs connecting the particles.
    dt (float):                     The time step for the simulation.
    uptime (float):                 The total time that has passed in the simulation.
    callback (Optional[Callable]):  An optional callback function for custom operations.
    """

    def __init__(self, dt=0.01, callback: Optional[Callable] = None):
        self.particles: List[Particle] = []
        self.springs: List[Spring] = []
        self.dt: float = dt
        self.uptime = 0.0

        self.callback = callback

    def apply_forces(self) -> None:
        """
        Applies forces to all particles in the simulation.

        This method first resets the forces for each particle, then applies
        spring forces, particle forces, and gravitational forces as necessary.
        """

        # 1. Reset Forces
        for particle in self.particles:
            particle.forces = []

        if self.callback:
            self.callback(self)

        # 2. Spring Forces
        for spring in self.springs:
            spring.apply_forces()

        # 3. Particle Forces
        for i in range(len(self.particles)):
            # apply Earth gravity forces
            if not self.particles[i].is_anchored:
                self.particles[i].apply_forces()

            # 4. Universal Gravitation Equations (if applicable)
            for j in range(len(self.particles)):
                if i < j:
                    apply_gravitational_forces(self.particles[i], self.particles[j])

    # time intergration methods
    def apply_symplectic_euler(self):
        """
        Updates the simulation state using the Symplectic Euler method.

        This method is a time integration technique that's particularly useful for physical simulations.
        It first updates the velocity of each particle based on its acceleration, and then updates
        the position based on the new velocity. It's a simple yet effective method for many simulations.
        """
        self.apply_forces()
        for particle in self.particles:
            particle.velocity += self.dt * compute_acceleration(particle)
            particle.position += self.dt * particle.velocity

    def apply_exlicit_euler(self):
        """
        Updates the simulation state using the Explicit Euler method.

        This method is a straightforward time integration technique. It first updates the position
        of each particle using its current velocity, and then updates the velocity based on the
        acceleration. This method is easy to understand and implement.
        """
        self.apply_forces()
        for particle in self.particles:
            particle.position += self.dt * particle.velocity
            particle.velocity += self.dt * compute_acceleration(particle)

    def apply_leapfrog(self):
        """
        Updates the simulation state using the Leapfrog method.

        The Leapfrog method is a time integration technique that is especially good for simulations
        requiring stable energy conservation over time. It updates velocities and positions in a staggered
        way, leading to improved numerical stability over some other methods.
        """

        # Step 1: Apply forces to compute accelerations
        self.apply_forces()

        # Step 2: Update velocities by half a timestep
        for particle in self.particles:
            if not particle.is_anchored:
                particle.velocity += 0.5 * self.dt * compute_acceleration(particle)

        # Step 3: Update positions
        for particle in self.particles:
            particle.position += self.dt * particle.velocity

        # Step 4: Apply forces again (as positions have changed)
        self.apply_forces()

        # Step 5: Complete the velocity update
        for particle in self.particles:
            if not particle.is_anchored:
                particle.velocity += 0.5 * self.dt * compute_acceleration(particle)

    def apply_runge_kutta_2(self):
        """
        Updates the simulation state using the second-order Runge-Kutta method (RK2).

        RK2 is a more sophisticated method for time integration, providing better accuracy
        than the simpler Euler methods. It updates the positions and velocities using an intermediate
        step to compute a more accurate final result.
        """

        initial_velocities = [particle.velocity for particle in self.particles]
        initial_positions = [particle.position for particle in self.particles]

        # Step 1: Apply forces to compute initial accelerations
        self.apply_forces()
        k1_velocities = [self.dt * compute_acceleration(particle) for particle in self.particles]
        k1_positions = [self.dt * particle.velocity for particle in self.particles]  # based on previous velocity

        # Step 3: Update forces based on midpoint positions
        for i, particle in enumerate(self.particles):
            if not particle.is_anchored:
                particle.velocity = initial_velocities[i] + 0.5 * k1_velocities[i]
                particle.position = initial_positions[i] + 0.5 * k1_positions[i]

        self.apply_forces()
        k2_velocities = [self.dt * compute_acceleration(particle) for particle in self.particles]
        k2_positions = [self.dt * (particle.velocity + 0.5 * k1_velocities[i]) for i, particle in
                        enumerate(self.particles)]

        # Step 4: Update velocity and position based on midpoint accelerations
        for i, particle in enumerate(self.particles):
            if not particle.is_anchored:
                particle.velocity = initial_velocities[i] + k2_velocities[i]
                particle.position = initial_positions[i] + k2_positions[i]

    def apply_runge_kutta_4(self):
        """
        Updates the simulation state using the fourth-order Runge-Kutta method (RK4).

        RK4 is one of the most commonly used methods for solving ordinary differential equations
        in simulations. It offers excellent accuracy through a series of intermediate steps, making it
        ideal for more complex or sensitive simulations where precision is key.
        """

        # save current state
        initial_velocities = [particle.velocity for particle in self.particles]  # added
        initial_positions = [particle.position for particle in self.particles]

        # Step 1: Apply forces and compute initial accelerations
        self.apply_forces()
        k1_velocities = [self.dt * compute_acceleration(particle) for particle in self.particles]
        k1_positions = [self.dt * particle.velocity for particle in self.particles]  # based on previous velocity

        # Step 2: Compute mid-point velocities and positions for k2
        for i, particle in enumerate(self.particles):
            if not particle.is_anchored:
                particle.velocity = initial_velocities[i] + 0.5 * k1_velocities[i]
                particle.position = initial_positions[i] + 0.5 * k1_positions[i]

        self.apply_forces()
        k2_velocities = [self.dt * compute_acceleration(particle) for particle in self.particles]
        k2_positions = [self.dt * (particle.velocity + 0.5 * k1_velocities[i]) for i, particle in
                        enumerate(self.particles)]

        # Step 3: Compute mid-point velocities and positions for k3
        for i, particle in enumerate(self.particles):
            if not particle.is_anchored:
                particle.velocity = initial_velocities[i] + 0.5 * k2_velocities[i]
                particle.position = initial_positions[i] + 0.5 * k2_positions[i]

        self.apply_forces()
        k3_velocities = [self.dt * compute_acceleration(particle) for particle in self.particles]
        k3_positions = [self.dt * (particle.velocity + 0.5 * k2_velocities[i]) for i, particle in
                        enumerate(self.particles)]

        # Step 4: Compute end velocities and positions for k4
        for i, particle in enumerate(self.particles):
            if not particle.is_anchored:
                particle.velocity = initial_velocities[i] + k3_velocities[i]
                particle.position = initial_positions[i] + k3_positions[i]

        self.apply_forces()
        k4_velocities = [self.dt * compute_acceleration(particle) for particle in self.particles]
        k4_positions = [self.dt * (particle.velocity + k3_velocities[i]) for i, particle in enumerate(self.particles)]

        # Step 5: Update final positions and velocities
        for i, particle in enumerate(self.particles):
            if not particle.is_anchored:
                particle.velocity = initial_velocities[i] + (
                        k1_velocities[i] + 2 * k2_velocities[i] + 2 * k3_velocities[i] + k4_velocities[i]) / 6
                particle.position = initial_positions[i] + (
                        k1_positions[i] + 2 * k2_positions[i] + 2 * k3_positions[i] + k4_positions[i]) / 6


# Simulation class
class Simulation:
    """
    Main class to control and run the physics simulation.

    This class encapsulates the simulation state and provides an interface to run the simulation
    with different integrators. It also handles callbacks for forces and statistics,
    and manages the visualization aspects like drawing and updating the objects.

    Attributes:
    state (SimulationState):            The current state of the simulation.
    callback (Callable, optional):      Function to call after each simulation step.
    stats_callback (Callable, optional): Function to call for generating statistics.
    integrator (str):                   The integrator method used for time integration.
    """

    def __init__(self, dt=0.1, callback=None, forces_callback=None, stats_callback=None,
                 integrator: str = INTEGRATOR_SYMPLECTIC_EULER):
        """
        Initializes the Simulation object.

        Parameters:
        dt (float):                             The time step for the simulation.
        callback (Callable, optional):          Function to call after each simulation step.
        forces_callback (Callable, optional):   Function to call for applying additional forces.
        stats_callback (Callable, optional):    Function to call for generating statistics.
        integrator (str):                  The integrator method used for time integration (default: symplectic-euler).
        """
        self.state = SimulationState(dt=dt, callback=forces_callback)
        self.callback = callback
        self.stats_callback = stats_callback
        self.integrator = integrator

        # internal profiling
        self.last_step_duration = 0.0
        self.last_time = time.time()

    @property
    def particles(self):
        return self.state.particles

    @property
    def springs(self):
        return self.state.springs

    @property
    def dt(self):
        return self.state.dt

    @property
    def uptime(self):
        return self.state.uptime

    @property
    def forces_callback(self):
        return self.state.callback

    @forces_callback.setter
    def forces_callback(self, callback):
        self.state.callback = callback

    def get_stats(self):
        """
        Gathers statistics about the current state of the simulation.

        This method compiles various statistical data about the simulation, like the current time,
        total energy, and the duration of the last simulation step. If a stats callback is set,
        it includes additional statistics from that as well.

        Returns:
        str: A formatted string containing relevant statistics about the simulation.
        """
        return "\n".join([
            f"$t={self.state.uptime:.2f}$",
            f"Energy={self.total_energy():.2f}",
            f"Step Duration={self.last_step_duration * 1000:.0f}ms",
        ]) + ("\n" + self.stats_callback(self) if self.stats_callback else "")

    def simulation_step(self):
        """
        Executes a single step of the simulation.

        This method updates the simulation state by one time step. Depending on the chosen
        integrator method, it applies the relevant physics and force updates. It also
        handles any defined callbacks and updates the simulation's uptime.
        """

        self.last_time = time.time()

        if self.callback:
            self.callback(self)

        # apply time integration
        if self.integrator == INTEGRATOR_SYMPLECTIC_EULER:
            self.state.apply_symplectic_euler()
        elif self.integrator == INTEGRATOR_EXPLICIT_EULER:
            self.state.apply_exlicit_euler()
        elif self.integrator == INTEGRATOR_LEAPFROG:
            self.state.apply_leapfrog()
        elif self.integrator == INTEGRATOR_RK2:
            self.state.apply_runge_kutta_2()
        elif self.integrator == INTEGRATOR_RK4:
            self.state.apply_runge_kutta_4()

        current_time = time.time()
        self.last_step_duration = current_time - self.last_time
        self.last_time = current_time

        self.state.uptime += self.state.dt

    #
    def get_artists(self) -> List[Artist]:
        """
        Retrieves the matplotlib artists needed for visualization.

        This method collects all the drawable objects (like particles and springs) that
        need to be displayed on the plot. It's particularly useful for updating the
        visualization in each frame of the animation.

        Returns:
        List[matplotlib.artist.Artist]: A list of drawable objects for the plot.
        """
        artists = []
        for spring in self.springs:
            artists += spring.get_artists()
        for particle in self.particles:
            artists += particle.get_artists()
        return artists

    def get_legend_handles(self) -> List[Artist]:
        """
        Retrieves legend handles for the plot.

        This method is used to create a legend for the plot. It goes through all the particles
        and collects their respective legend handles if they exist.

        Returns:
        List[matplotlib.artist.Artist]: A list of legend handles for the particles.
        """
        artists = []
        artists += [particle.get_legend_handle() for particle in self.particles if particle.get_legend_handle()]
        return artists

    def update_artists(self):
        """
        Updates the position and visuals of the drawable objects.

        This method is called to refresh the positions and states of all the drawable objects
        (particles and springs) in the plot, ensuring the visualization reflects the current
        state of the simulation.
        """

        for particle in self.particles:
            particle.update_artist()
        for spring in self.springs:
            spring.update_artist()

    def get_force_quivers(self) -> Tuple[List, List, List]:
        """
        Retrieves the data for drawing force vectors as quivers.

        This method is used to visualize the forces acting on each particle. It compiles
        the data necessary to draw quivers (arrows) representing these forces.

        Returns:
        Tuple[List, List, List]: A tuple containing lists of offsets, horizontal components,
                                 and vertical components of the quiver arrows.
        """
        offsets, us, vs = [], [], []
        for particle in self.particles:
            p_offsets, p_us, p_vs = particle.get_forces_quivers()
            offsets += p_offsets
            us += p_us
            vs += p_vs
        return offsets, us, vs

    def total_energy(self) -> float:
        """
        Calculates the total energy of the system.

        This method computes the total energy of the system by summing up the kinetic energy
        of the particles, the potential energy stored in springs, and the gravitational
        potential energy (if applicable).

        Returns:
        float: The total energy of the system.
        """

        total_ke = 0.0          # Total kinetic energy
        total_pe_spring = 0.0   # Total potential energy in springs
        total_pe_gravity = 0.0  # Total gravitational potential energy

        # Calculate total kinetic energy
        for particle in self.particles:
            if not particle.is_anchored:
                total_ke += particle.get_kinetic_energy()

        # Calculate total potential energy in springs
        for spring in self.springs:
            total_pe_spring += spring.get_potential_energy()

        total_energy = total_ke + total_pe_spring + total_pe_gravity
        return total_energy
