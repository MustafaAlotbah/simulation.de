import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches
from typing import *


earth_gravity = np.array([0, -9.81], dtype='float64')
gravitational_constant = 6.67384 * (10 ** -6)  # This constant is NOT realistic!

class Particle:
    def __init__(self,
                 x=0,
                 y=0,
                 vx=0,
                 vy=0,
                 mass=1,
                 color='white',
                 label="",
                 has_legend=True,
                 is_anchored=False,
                 trace_alpha=0.6,
                 radius=None,
                 display_forces=True
                 ):
        """
        Represents a physical particle in the simulation.

        The particle has physical properties like position, velocity, and mass, along with
        visual properties for simulation visualization, such as color and radius.

        Parameters:
        x (float):                  The initial x-coordinate of the particle.
        y (float):                  The initial y-coordinate of the particle.
        vx (float):                 The initial velocity of the particle along the x-axis.
        vy (float):                 The initial velocity of the particle along the y-axis.
        mass (float):               The mass of the particle.
        color (str):                The color of the particle in the visualization.
        label (str):                Label for the particle, used in the legend.
        has_legend (bool):          Indicates if the particle should appear in the legend.
        is_anchored (bool):         If True, the particle does not move.
        trace_alpha (float):        The transparency of the trace line.
        radius (float, optional):   The radius of the particle in the visualization.
        display_forces (bool):      If True, displays force vectors acting on the particle.
        """

        # physical properties
        self.position = np.array([x, y], dtype=np.float)
        self.velocity = np.array([vx, vy], dtype=np.float)
        self.forces = []
        self.mass = mass

        # visuals
        self.alpha = 1.0
        self.trace_alpha = trace_alpha
        self.has_legend = has_legend
        self.display_forces = display_forces

        # default radius is handcrafted and is not important
        if not radius:
            radius = max(1, np.sqrt(max(0, np.log10(self.mass)))) / 30
        self.radius = radius

        # other properties
        self.color = color
        self.label = label
        self.is_anchored = is_anchored
        self.apply_earth_gravity = True
        self.apply_gravitational_forces = False  # This is only for the simulation of orbits

        # what to draw on the plot
        self.artist = plt.Circle(
            tuple(self.position),
            radius=1,
            fc=self.color, ec='black',
            linewidth=2,
            zorder=10
        )

        # what to draw on the legend
        self.legend_handle = plt.Line2D(  # this is a -o- plot
            [0], [0],
            color='#ffffff00',  # this is the line's color (transparent)
            markerfacecolor=self.color,
            markeredgecolor=self.color,
            marker='o',
            markersize=10,
            label=self.label
        )

        # information to keep track of its position
        self.positions = [tuple(self.position)]

        _path = mpl_patches.Path(self.positions, [mpl_patches.Path.MOVETO])
        self.path_artist = mpl_patches.PathPatch(
            path=_path, ec=self.color, lw=0.75, fill=False, zorder=8, alpha=self.trace_alpha
        )

    def get_legend_handle(self) -> Optional[plt.Line2D]:
        """
        Retrieves the legend handle for the particle.

        Returns the handle used to represent this particle in the plot's legend, if applicable.

        Returns:
        matplotlib.patches.Line2D or None: The legend handle or None if the particle has no legend.
        """
        if self.has_legend:
            return self.legend_handle
        return None

    def get_artists(self) -> Tuple[plt.Circle, mpl_patches.PathPatch]:
        """
        Retrieves the matplotlib artists for the particle.

        These artists are used to draw the particle and its trace in the plot.

        Returns:
        Tuple[matplotlib.patches.Circle, matplotlib.patches.PathPatch]:
            The circle representing the particle and the path representing its trace.
        """
        return self.artist, self.path_artist

    def update_artist(self) -> None:
        """
        Updates the particle's artists based on the current state.

        This method is called to update the visual representation of the particle,
        including its position and trace.
        """

        # update the particle itself
        self.artist.center = self.position
        self.artist.set_facecolor(self.color)
        self.artist.radius = self.radius  # size by design
        self.artist.set_alpha(self.alpha)

        # update the particle's trace
        self.positions.append(tuple(self.position))
        self.positions = self.positions[-750:]  # keep the tail 750 points long at most
        _path = mpl_patches.Path(
            self.positions, [mpl_patches.Path.MOVETO] + [mpl_patches.Path.LINETO] * (len(self.positions) - 1)
        )
        self.path_artist.set_path(_path)

    def get_forces_quivers(self) -> Tuple[List[Tuple[float, float]], List[float], List[float]]:
        """
        Retrieves the data for visualizing the forces acting on the particle.

        This method compiles the data necessary for drawing quivers (arrows) representing the forces.

        Returns:
        Tuple[List, List, List]: Lists of offsets, horizontal components, and vertical components for the force quivers.
        """

        # maybe the set of force
        if not self.display_forces:
            return [], [], []
        forces_offset: List[Tuple[float, float]] = [self.position.tolist() for _ in self.forces]
        forces_u = [force[0] for force in self.forces]
        forces_v = [force[1] for force in self.forces]
        return forces_offset, forces_u, forces_v

    def apply_forces(self) -> None:
        """
        Applies the relevant forces to the particle.

        This method updates the force vector of the particle based on the simulation conditions,
        such as applying gravity.
        """

        if self.apply_earth_gravity:
            self.forces.append(earth_gravity * self.mass)

    def get_kinetic_energy(self) -> float:
        """
        Calculates and returns the kinetic energy of the particle.

        Returns:
        float: The kinetic energy of the particle.
        """
        speed = np.linalg.norm(self.velocity)
        kinetic_energy: float = 0.5 * self.mass * speed ** 2
        gravity_kinetic_energy: float = 0

        if self.apply_earth_gravity:
            gravity_kinetic_energy = self.mass * 9.81 * self.position[1]  # y axis

        return kinetic_energy + gravity_kinetic_energy


class Spring:
    def __init__(self, particle1: Particle, particle2: Particle, stiffness: Optional[float] = None, length=None, damping=1, dt=0.005):
        """
        Represents a spring connecting two particles in the simulation.

        The spring exerts forces on the particles based on its length and stiffness,
        simulating elastic connections in the physical system.

        Parameters:
        particle1 (Particle):   The first particle connected by the spring.
        particle2 (Particle):   The second particle connected by the spring.
        stiffness (float):      The stiffness of the spring.
        length (float, optional): The rest length of the spring. If None, the initial distance between the particles is used.
        damping (float):        The damping coefficient of the spring.
        dt (float):             The time step for the simulation, used in stiffness calculation if stiffness is None.
        """

        self.particle1 = particle1
        self.particle2 = particle2

        # default length is the distance between both particles when the spring was created
        if not length:
            length = np.linalg.norm(particle1.position - particle2.position)

        if not stiffness:
            stiffness = 1/(
                    dt**2 * (1/particle1.mass + 1/particle2.mass)
            )

        self.rest_length = length
        self.stiffness = stiffness
        self.damping = damping

        # what to draw on the plot
        self._path = mpl_patches.Path(
            [self.particle1.position, self.particle2.position],
            [mpl_patches.Path.MOVETO, mpl_patches.Path.LINETO]
        )
        self.artist = mpl_patches.PathPatch(path=self._path, ec='white', lw=1.75, fill=False, zorder=8, alpha=1.0)

    def get_artists(self) -> Tuple[mpl_patches.PathPatch,]:
        """
        Retrieves the matplotlib artist for the spring.

        This artist is used to draw the spring in the plot.

        Returns:
        Tuple[matplotlib.patches.PathPatch]: The path patch representing the spring.
        """
        return self.artist,

    def update_artist(self) -> None:
        """
        Updates the spring's artist based on the current state.

        This method is called to update the visual representation of the spring,
        adjusting its position to match the connected particles.
        """
        # update the spring itself
        self._path.vertices = np.array([self.particle1.position, self.particle2.position])

    def get_potential_energy(self) -> float:
        """
        Calculates and returns the potential energy stored in the spring.

        Returns:
        float: The potential energy of the spring based on its extension or compression.
        """
        current_length = np.linalg.norm(self.particle1.position - self.particle2.position)
        spring_extension = current_length - self.rest_length
        potential_energy = 0.5 * self.stiffness * spring_extension ** 2
        return potential_energy

    def apply_forces(self) -> None:
        """
        Applies the spring forces to the connected particles.

        This method calculates the elastic and damping forces based on the spring's
        current length, rest length, and the velocities of the connected particles.
        """

        # current length of the spring as a vector
        x = self.particle1.position - self.particle2.position

        # current length of the spring as a scalar
        norm_x = np.linalg.norm(x)

        # compute the displacement: ||x1 - x2|| - restLength
        C = norm_x - self.rest_length

        if norm_x > 0.001:
            # gradient of the Constraint C:
            # dC/dX = (X1 - X2) / || X1 - X2) ||
            grad_c = x / norm_x

            # - Spring force :
            # F_spring = - k . C . dC/dX
            F_spring = -self.stiffness * C * grad_c

            # \dot{C} = (V1 - V2) / (dC/dX)
            gradc_wrt_t = np.dot(self.particle1.velocity - self.particle2.velocity, grad_c)

            # - Damping force
            # F_damping = - k' . \dot{C} . dC/dX
            F_damping = -self.damping * gradc_wrt_t * grad_c

            if not self.particle1.is_anchored:
                self.particle1.forces += [F_spring + F_damping]

            if not self.particle2.is_anchored:
                self.particle2.forces += [-(F_spring + F_damping)]
