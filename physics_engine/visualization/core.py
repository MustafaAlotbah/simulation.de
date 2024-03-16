"""
    Core of the visualization module for simulation.de visualizations.
    This module provides the functionality for visualizing the simulations created with the physics engine.
    Author: Mustafa Alotbah, 2024
"""

import numpy as np
import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.pyplot as plt

from physics_engine.simulation import Simulation

import time


class Visualization:
    """
    Initializes the visualization for the simulation.

    This class sets up the matplotlib figure and axes for plotting, applies styling, and
    initializes the animation for visualizing the simulation steps.

    Parameters:
    simulation (Simulation):    The simulation object to visualize.
    bounds (list, optional):    The bounds of the visualization area as [xmin, xmax, ymin, ymax].
    simulation_steps (int):     The number of simulation steps to visualize.
    foreground_color (str, optional): The color of the foreground elements (text, ticks).
    background_color (str, optional): The background color of the plot.
    frames_skipped (int):       The number of frames to skip between updates in the visualization.
    """
    def __init__(self,
                 simulation: Simulation,
                 bounds=None,
                 simulation_steps=2501,
                 foreground_color=None,
                 background_color=None,
                 frames_skipped=0
                 ):
        assert frames_skipped >= 0  # Ensures that the number of frames skipped is non-negative
        self.frames_skipped = frames_skipped

        # Set default foreground and background colors if not provided
        self.foreground_color = foreground_color if foreground_color else mpl.rcParams["text.color"]
        self.background_color = background_color if background_color else mpl.rcParams["axes.labelcolor"]

        # Set the grid color for the plot
        self.grid_color = '#181818'

        # Apply color settings to matplotlib's global parameters
        mpl.rcParams["text.color"] = self.foreground_color
        mpl.rcParams["axes.labelcolor"] = self.background_color
        mpl.rcParams["xtick.color"] = self.foreground_color
        mpl.rcParams["ytick.color"] = self.foreground_color

        # Set the bounds for the plot area
        self.bounds = bounds if bounds else [-5, 5, -10, 10]

        # Reference to the simulation object to be visualized
        self.simulation = simulation

        # Initialize the matplotlib figure and axes for the visualization
        self.fig, self.ax1 = plt.subplots(1, 1, tight_layout=True)
        self.fig.set_facecolor(self.background_color)
        self.fig.set_size_inches(9, 16, forward=True)

        # Configure the axes, including titles and colors
        self.ax1.set_title("Space")
        self.ax1.set_facecolor(self.background_color)

        # Add a legend to the plot with simulation entities
        self.ax1.legend(
            handles=self.simulation.get_legend_handles(),
            framealpha=0.3,
            bbox_to_anchor=(0.98, 0.93), loc='upper right'
        )

        # Set up the plot boundaries and aspect ratio
        self.ax1.axis(self.bounds)
        self.ax1.set_aspect(1)

        # Configure the plot's ticks and grid
        self.ax1.set_xticks(range(self.bounds[0], self.bounds[1]))
        self.ax1.set_yticks(range(self.bounds[2], self.bounds[3]))
        self.ax1.grid(c=self.grid_color, zorder=-10)

        # Add a text object to the plot for displaying statistics
        self.stats_text = self.ax1.text(0.02, 0.92, "stats",
                                        zorder=100,
                                        transform=self.ax1.transAxes,
                                        verticalalignment='top'
                                        )
        self.stats_text.set_bbox(dict(facecolor='0.80', alpha=0.3, edgecolor=background_color))

        # Initialize a quiver plot for visualizing forces
        self.forces_quiver = self.ax1.quiver([0], [0], [1], [1], cmap="winter", width=0.005, alpha=0.9, zorder=100)

        # Add artists for each entity in the simulation
        for artist in self.simulation.get_artists():
            self.ax1.add_patch(artist)

        def anim_update(_):
            """ Function to update the animation at each step """
            start_time = time.time()

            # Update artists and perform simulation steps
            simulation.update_artists()
            for i in range(frames_skipped+1):
                simulation.simulation_step()

            # Update the statistics text in the plot
            self.stats_text.set_text(simulation.get_stats())

            # Get and normalize the force vectors for quivers
            offsets, us, vs = simulation.get_force_quivers()
            assert len(offsets) == len(us) and len(offsets) == len(vs)

            us = np.array(us)
            vs = np.array(vs)
            magnitudes = np.sqrt(us ** 2 + vs ** 2)
            normalized_us = us / magnitudes * 2
            normalized_vs = vs / magnitudes * 2
            thicknesses = np.log(magnitudes + 1)

            # Update the quiver plot with the new force vectors
            self.forces_quiver.N = len(offsets)
            self.forces_quiver.set_offsets(offsets)
            self.forces_quiver.set_UVC(normalized_us, normalized_vs, thicknesses)

            # Print the frame number and processing time periodically
            if _ % 10 == 0:
                print(_, f"({self.simulation.uptime:.3f} s)", f"~ {(time.time() - start_time)*1000:.0f} ms")

            return self.stats_text, self.forces_quiver, *simulation.get_artists()

        self.animation = mpl.animation.FuncAnimation(
            fig=self.fig,
            func=anim_update,
            interval=10,
            blit=True,
            frames=simulation_steps,
            repeat=False
        )

    @staticmethod
    def show():
        """
        Displays the matplotlib plot.

        This method should be called to actually display the plot window after setting up
        the visualization.
        """
        plt.show()

    def save_video(self, file_name, dpi=300):
        """
        Saves the visualization as a high-quality video file.

        The video resolution is determined by the DPI (dots per inch) parameter. Higher DPI values
        result in higher resolution videos.

        Parameters:
        file_name (str):    The path and file name where the video will be saved.
        dpi (int, optional): Dots per inch for the video resolution. Default is 300 for high quality.
        """
        fps = 1/(self.simulation.dt * (1+self.frames_skipped))
        print(f"FPS: {fps:.2f}")
        self.animation.save(
            file_name,
            fps=fps,
            dpi=dpi,
            extra_args=['-vcodec', 'libx264'],
            savefig_kwargs={'facecolor': self.background_color}
        )
